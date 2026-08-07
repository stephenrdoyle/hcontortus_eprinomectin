#!/usr/bin/env bash
#
# PREREQUISITES
# - you are on an HPC cluster
# - your job scheduler is LSF (bsub)
# - Singularity module is loaded
# - the following module is installed and available by module load : nextflow
# - the following module is installed and available by module load : mapping-helminth/v1.0.8
# - in the current path, you created a folder 0_refseq, in which you donwloaded refseq files .fa (reference genome) and .gtf (gene annotation)
# - the .manifest file is present in the current launch directory, and looks like :
#   ID,R1,R2
#   sample_1,ftp://.../sample_1_R1.fastq.gz,ftp://.../sample_1_R2.fastq.gz
#
# EXECUTION PRINCIPLE
# - if an analysis-stage directory is absent, a temporary stage directory is created and launched,
# - the temporary stage directory is renamed to the final stage directory only after successful completion,
# - if the run fails, manually remove the remaining temporary stage directory and rerun,
# - each output is first written to a .tmp.XXXXXX path and is renamed only after successful completion, so remaining .tmp.XXXXXX files or directories mean that the stage is corrupted.
#
# WARNING
# Launch this pipeline using bsub on a long queue to avoid any interruption before it’s terminated.
# Command example : bsub -q long -J genomics -n 1 -R "select[mem>=2000] rusage[mem=2000] span[hosts=1]" -M 2000 -o "genomics.out" -e "genomics.err" bash genomics.sh

set -euo pipefail
IFS=$'\n\t'

ROOT="$(pwd -P)"
MANIFEST="$(find "${ROOT}" -maxdepth 1 -type f -name '*.manifest' -print -quit)"
REFERENCE="$(find "${ROOT}/0_refseq" -maxdepth 1 -type f -name '*.fa' -print -quit)"
ANNOTATION="$(find "${ROOT}/0_refseq" -maxdepth 1 -type f -name '*.gtf' -print -quit)"

# Wait for an LSF job and stop the pipeline if the job did not finish
# successfully. After bwait reports that the job has ended, wait until
# bjobs exposes its final DONE or EXIT status.
wait_for_job() {
    local jid="$1"
    local status

    if [[ ! "${jid}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: invalid or missing LSF job ID: '${jid}'." >&2
        exit 1
    fi

    bwait -w "ended(${jid})"

    while true; do
        status="$(
            bjobs -a -noheader -o "stat" "${jid}" 2>/dev/null \
            | awk 'NF {print $1; exit}' \
            || true
        )"

        case "${status}" in
            DONE)
                return 0
                ;;
            EXIT)
                echo "ERROR: LSF job ${jid} ended with status EXIT." >&2
                echo "Pipeline stopped. Remove the failed temporary stage directory before rerunning." >&2
                exit 1
                ;;
            *)
                sleep 2
                ;;
        esac
    done
}

# ==============================================================================
# 0_samples — sequential acquisition of paired-end FASTQ files with wget -c
# ==============================================================================
if [[ ! -d "${ROOT}/0_samples" ]]; then
    echo "Creating and running 0_samples"

    STAGE="$(mktemp -d "${ROOT}/0_samples.tmp.XXXXXX")"

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        r1="${r1//$'\r'/}"
        r2="${r2//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue

        mkdir -p "${STAGE}/${id}"

        submission="$(
            ROOT="${ROOT}" STAGE="${STAGE}" ID="${id}" URL_R1="${r1}" URL_R2="${r2}" \
            bsub \
                -J "DL_${id}" \
                -n 1 \
                -R "select[mem>=2000] rusage[mem=2000] span[hosts=1]" \
                -M 2000 \
                -o "${STAGE}/${id}/DL_${id}.%J.out" \
                -e "${STAGE}/${id}/DL_${id}.%J.err" \
                bash -lc '
set -euo pipefail

final_r1="${STAGE}/${ID}/${ID}_R1.fastq.gz"
tmp_r1="$(find "${STAGE}/${ID}" -maxdepth 1 -type f -name ".${ID}_R1.fastq.gz.tmp.*" -print -quit)"
if [[ -z "${tmp_r1}" ]]; then
    tmp_r1="$(mktemp "${STAGE}/${ID}/.${ID}_R1.fastq.gz.tmp.XXXXXX")"
fi
wget -c -O "${tmp_r1}" "${URL_R1}"
mv "${tmp_r1}" "${final_r1}"

final_r2="${STAGE}/${ID}/${ID}_R2.fastq.gz"
tmp_r2="$(find "${STAGE}/${ID}" -maxdepth 1 -type f -name ".${ID}_R2.fastq.gz.tmp.*" -print -quit)"
if [[ -z "${tmp_r2}" ]]; then
    tmp_r2="$(mktemp "${STAGE}/${ID}/.${ID}_R2.fastq.gz.tmp.XXXXXX")"
fi
wget -c -O "${tmp_r2}" "${URL_R2}"
mv "${tmp_r2}" "${final_r2}"
'
        )"

        echo "${submission}"
        jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
        wait_for_job "${jid}"
    done < "${MANIFEST}"

    mv "${STAGE}" "${ROOT}/0_samples"
fi


# ==============================================================================
# 1_QC — parallel FastQC assessment followed by MultiQC summarisation
# ==============================================================================
if [[ ! -d "${ROOT}/1_QC" ]]; then
    echo "Creating and running 1_QC"

    STAGE="$(mktemp -d "${ROOT}/1_QC.tmp.XXXXXX")"

    mkdir -p "${ROOT}/0_tools"

    if [[ ! -f "${ROOT}/0_tools/fastqc_0.12.1.sif" ]]; then
        fastqc_tmp="$(find "${ROOT}/0_tools" -maxdepth 1 -type f -name '.fastqc_0.12.1.sif.tmp.*' -print -quit)"
        if [[ -z "${fastqc_tmp}" ]]; then
            fastqc_tmp="$(mktemp "${ROOT}/0_tools/.fastqc_0.12.1.sif.tmp.XXXXXX")"
        fi
        wget -c \
            -O "${fastqc_tmp}" \
            "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
        mv "${fastqc_tmp}" "${ROOT}/0_tools/fastqc_0.12.1.sif"
    fi

    if [[ ! -f "${ROOT}/0_tools/multiqc_1.17.sif" ]]; then
        multiqc_tmp="$(find "${ROOT}/0_tools" -maxdepth 1 -type f -name '.multiqc_1.17.sif.tmp.*' -print -quit)"
        if [[ -z "${multiqc_tmp}" ]]; then
            multiqc_tmp="$(mktemp "${ROOT}/0_tools/.multiqc_1.17.sif.tmp.XXXXXX")"
        fi
        wget -c \
            -O "${multiqc_tmp}" \
            "https://depot.galaxyproject.org/singularity/multiqc%3A1.17--pyhdfd78af_1"
        mv "${multiqc_tmp}" "${ROOT}/0_tools/multiqc_1.17.sif"
    fi

    mkdir -p "${STAGE}/fastqc"
    jobs=()

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue

        submission="$(
            ROOT="${ROOT}" STAGE="${STAGE}" ID="${id}" \
            bsub \
                -J "QC_${id}" \
                -n 4 \
                -R "select[mem>=4000] rusage[mem=4000] span[hosts=1]" \
                -M 4000 \
                -o "${STAGE}/fastqc/QC_${id}.%J.out" \
                -e "${STAGE}/fastqc/QC_${id}.%J.err" \
                bash -lc '
set -euo pipefail

tmpdir="$(mktemp -d "${STAGE}/fastqc/.${ID}.tmp.XXXXXX")"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/fastqc_0.12.1.sif" \
    fastqc \
    --threads 4 \
    --outdir "${tmpdir}" \
    "${ROOT}/0_samples/${ID}/${ID}_R1.fastq.gz" \
    "${ROOT}/0_samples/${ID}/${ID}_R2.fastq.gz"

mv "${tmpdir}" "${STAGE}/fastqc/${ID}"
'
        )"

        echo "${submission}"
        jobs+=("$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")")
    done < "${MANIFEST}"

    for jid in "${jobs[@]}"; do
        wait_for_job "${jid}"
    done

    submission="$(
        ROOT="${ROOT}" STAGE="${STAGE}" \
        bsub \
            -J "QC_MultiQC" \
            -n 1 \
            -R "select[mem>=8000] rusage[mem=8000] span[hosts=1]" \
            -M 8000 \
            -o "${STAGE}/QC_MultiQC.%J.out" \
            -e "${STAGE}/QC_MultiQC.%J.err" \
            bash -lc '
set -euo pipefail

tmpdir="$(mktemp -d "${STAGE}/.multiqc.tmp.XXXXXX")"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/multiqc_1.17.sif" \
    multiqc \
    --filename raw_reads_multiqc.html \
    --outdir "${tmpdir}" \
    "${STAGE}/fastqc"

mv "${tmpdir}" "${STAGE}/multiqc"
'
    )"

    echo "${submission}"
    jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
    wait_for_job "${jid}"

    mv "${STAGE}" "${ROOT}/1_QC"
fi


# ==============================================================================
# 2_mapping — mapping-helminth v1.0.8 followed by MultiQC summarisation of
# SAMtools flagstat and stats outputs generated by the mapping pipeline
# ==============================================================================
if [[ ! -d "${ROOT}/2_mapping" ]]; then
    echo "Creating and running 2_mapping"

    STAGE="$(mktemp -d "${ROOT}/2_mapping.tmp.XXXXXX")"

    mkdir -p "${ROOT}/0_tools"

    if [[ ! -f "${ROOT}/0_tools/multiqc_1.17.sif" ]]; then
        multiqc_tmp="$(find "${ROOT}/0_tools" -maxdepth 1 -type f -name '.multiqc_1.17.sif.tmp.*' -print -quit)"
        if [[ -z "${multiqc_tmp}" ]]; then
            multiqc_tmp="$(mktemp "${ROOT}/0_tools/.multiqc_1.17.sif.tmp.XXXXXX")"
        fi
        wget -c \
            -O "${multiqc_tmp}" \
            "https://depot.galaxyproject.org/singularity/multiqc%3A1.17--pyhdfd78af_1"
        mv "${multiqc_tmp}" "${ROOT}/0_tools/multiqc_1.17.sif"
    fi

    manifest_tmp="$(mktemp "${STAGE}/.fqpaths.manifest.tmp.XXXXXX")"
    printf 'ID,R1,R2\n' > "${manifest_tmp}"

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue

        printf '%s,%s,%s\n' \
            "${id}" \
            "${ROOT}/0_samples/${id}/${id}_R1.fastq.gz" \
            "${ROOT}/0_samples/${id}/${id}_R2.fastq.gz" \
            >> "${manifest_tmp}"
    done < "${MANIFEST}"

    mv "${manifest_tmp}" "${STAGE}/fqpaths.manifest"

    submission="$(
        ROOT="${ROOT}" STAGE="${STAGE}" REFERENCE="${REFERENCE}" \
        bsub \
            -J "HcoMap" \
            -n 1 \
            -R "select[mem>=80000] rusage[mem=80000] span[hosts=1]" \
            -M 80000 \
            -o "${STAGE}/HcoMap.%J.out" \
            -e "${STAGE}/HcoMap.%J.err" \
            bash -lc '
set -euo pipefail

source /etc/profile.d/modules.sh
module load mapping-helminth/v1.0.8

# The mapping-helminth module currently exposes an obsolete shared Singularity
# image library (/software/pathogen/images). Disable this library and use a
# persistent project-local cache on the shared filesystem instead.
unset NEXTFLOW_SINGULARITY_LIBRARY
unset NXF_SINGULARITY_LIBRARYDIR

export NXF_SINGULARITY_CACHEDIR="${ROOT}/0_tools/nextflow_singularity_cache"
mkdir -p "${NXF_SINGULARITY_CACHEDIR}"

tmpdir="$(mktemp -d "${STAGE}/.outputs.tmp.XXXXXX")"
export NXF_WORK="${STAGE}/work"
export NXF_LOG_FILE="${STAGE}/nextflow.log"
mkdir -p "${NXF_WORK}" "${STAGE}/index_cache"

# mapping-helminth currently ignores some failed Nextflow processes.
# Force the complete Nextflow run to return a non-zero exit status whenever
# at least one process failed and was ignored.
cat > "${STAGE}/nextflow_override.config" <<EOF
workflow.failOnIgnore = true

singularity {
    cacheDir = "${ROOT}/0_tools/nextflow_singularity_cache"
}
EOF

mapping-helminth \
    -c "${STAGE}/nextflow_override.config" \
    --reference "${REFERENCE}" \
    --input "${STAGE}/fqpaths.manifest" \
    --outdir "${tmpdir}" \
    --index_cache "${STAGE}/index_cache" \
    --keep_unmapped true

# Independently verify that mapping produced the expected BAM and BAM index
# for every sample before accepting the stage as successfully completed.
missing=0

while IFS=, read -r id r1 r2; do
    id="${id//$'\r'/}"
    [[ "${id}" == "ID" || -z "${id}" ]] && continue

    bam="${tmpdir}/${id}/${id}.bam"

    if [[ ! -s "${bam}" ]]; then
        echo "ERROR: missing BAM for sample ${id}: ${bam}" >&2
        missing=1
    fi

    if [[ ! -s "${bam}.bai" && ! -s "${tmpdir}/${id}/${id}.bai" ]]; then
        echo "ERROR: missing BAM index for sample ${id}." >&2
        missing=1
    fi
done < "${STAGE}/fqpaths.manifest"

if (( missing != 0 )); then
    echo "ERROR: mapping completed incompletely. Pipeline stopped." >&2
    exit 1
fi

mv "${tmpdir}" "${STAGE}/outputs"
'
    )"

    echo "${submission}"
    jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
    wait_for_job "${jid}"

    submission="$(
        ROOT="${ROOT}" STAGE="${STAGE}" \
        bsub \
            -J "Map_MultiQC" \
            -n 1 \
            -R "select[mem>=8000] rusage[mem=8000] span[hosts=1]" \
            -M 8000 \
            -o "${STAGE}/Map_MultiQC.%J.out" \
            -e "${STAGE}/Map_MultiQC.%J.err" \
            bash -lc '
set -euo pipefail

tmpdir="$(mktemp -d "${STAGE}/.multiqc.tmp.XXXXXX")"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/multiqc_1.17.sif" \
    multiqc \
    --filename mapping_multiqc.html \
    --outdir "${tmpdir}" \
    "${STAGE}/outputs"

mv "${tmpdir}" "${STAGE}/multiqc"
'
    )"

    echo "${submission}"
    jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
    wait_for_job "${jid}"

    mv "${STAGE}" "${ROOT}/2_mapping"
fi


# ==============================================================================
# 3_variants — strict MAPQ >20 filtering, BCFtools v1.23 pile-up generation,
# multiallelic calling, biallelic normalisation and SnpEff v5.1d annotation
# SNPs and short insertions/deletions are retained.
# ==============================================================================
if [[ ! -d "${ROOT}/3_variants" ]]; then
    echo "Creating and running 3_variants"

    STAGE="$(mktemp -d "${ROOT}/3_variants.tmp.XXXXXX")"

    mkdir -p "${ROOT}/0_tools"

    if [[ ! -f "${ROOT}/0_tools/samtools_1.22.sif" ]]; then
        samtools_tmp="$(find "${ROOT}/0_tools" -maxdepth 1 -type f -name '.samtools_1.22.sif.tmp.*' -print -quit)"
        if [[ -z "${samtools_tmp}" ]]; then
            samtools_tmp="$(mktemp "${ROOT}/0_tools/.samtools_1.22.sif.tmp.XXXXXX")"
        fi
        wget -c \
            -O "${samtools_tmp}" \
            "https://depot.galaxyproject.org/singularity/samtools%3A1.22--h96c455f_0"
        mv "${samtools_tmp}" "${ROOT}/0_tools/samtools_1.22.sif"
    fi

    if [[ ! -f "${ROOT}/0_tools/bcftools_1.23.sif" ]]; then
        bcftools_tmp="$(find "${ROOT}/0_tools" -maxdepth 1 -type f -name '.bcftools_1.23.sif.tmp.*' -print -quit)"
        if [[ -z "${bcftools_tmp}" ]]; then
            bcftools_tmp="$(mktemp "${ROOT}/0_tools/.bcftools_1.23.sif.tmp.XXXXXX")"
        fi
        wget -c \
            -O "${bcftools_tmp}" \
            "https://depot.galaxyproject.org/singularity/bcftools%3A1.23--h3a4d415_0"
        mv "${bcftools_tmp}" "${ROOT}/0_tools/bcftools_1.23.sif"
    fi

    if [[ ! -f "${ROOT}/0_tools/snpeff_5.1d.sif" ]]; then
        snpeff_tmp="$(find "${ROOT}/0_tools" -maxdepth 1 -type f -name '.snpeff_5.1d.sif.tmp.*' -print -quit)"
        if [[ -z "${snpeff_tmp}" ]]; then
            snpeff_tmp="$(mktemp "${ROOT}/0_tools/.snpeff_5.1d.sif.tmp.XXXXXX")"
        fi
        wget -c \
            -O "${snpeff_tmp}" \
            "https://depot.galaxyproject.org/singularity/snpeff%3A5.1d--hdfd78af_2"
        mv "${snpeff_tmp}" "${ROOT}/0_tools/snpeff_5.1d.sif"
    fi

    if [[ ! -f "${REFERENCE}.fai" ]]; then
        reference_tmp="$(mktemp "$(dirname "${REFERENCE}")/.$(basename "${REFERENCE}").tmp.XXXXXX")"
        rm -f "${reference_tmp}"
        ln "${REFERENCE}" "${reference_tmp}"

        singularity exec \
            --bind "${ROOT}:${ROOT}" \
            "${ROOT}/0_tools/samtools_1.22.sif" \
            samtools faidx \
            "${reference_tmp}"

        rm "${reference_tmp}"
        mv "${reference_tmp}.fai" "${REFERENCE}.fai"
    fi

    mkdir -p "${STAGE}/filtered_bam"
    mkdir -p "${STAGE}/per_sample"

    # The species-specific database is built from the same WBPS18 genome and
    # annotation used for mapping and subsequent genomic interpretation.
    if [[ ! -d "${ROOT}/0_tools/snpeff_Hco_WBPS18" ]]; then
        submission="$(
            ROOT="${ROOT}" REFERENCE="${REFERENCE}" ANNOTATION="${ANNOTATION}" \
            bsub \
                -J "HcoSnpEffDB" \
                -n 4 \
                -R "select[mem>=12000] rusage[mem=12000] span[hosts=1]" \
                -M 12000 \
                -o "${ROOT}/0_tools/HcoSnpEffDB.%J.out" \
                -e "${ROOT}/0_tools/HcoSnpEffDB.%J.err" \
                bash -lc '
set -euo pipefail

tmpdir="$(mktemp -d "${ROOT}/0_tools/.snpeff_Hco_WBPS18.tmp.XXXXXX")"
mkdir -p "${tmpdir}/data/Hco_WBPS18"

cp "${REFERENCE}" \
   "${tmpdir}/data/Hco_WBPS18/sequences.fa"

cp "${ANNOTATION}" \
   "${tmpdir}/data/Hco_WBPS18/genes.gtf"

printf "data.dir = %s\nHco_WBPS18.genome : Haemonchus contortus PRJEB506 WBPS18\n" \
    "${tmpdir}/data" \
    > "${tmpdir}/snpEff.build.config"

printf "data.dir = %s\nHco_WBPS18.genome : Haemonchus contortus PRJEB506 WBPS18\n" \
    "${ROOT}/0_tools/snpeff_Hco_WBPS18/data" \
    > "${tmpdir}/snpEff.config"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/snpeff_5.1d.sif" \
    snpEff build \
    -gtf22 \
    -noCheckCds \
    -noCheckProtein \
    -c "${tmpdir}/snpEff.build.config" \
    -v Hco_WBPS18

rm "${tmpdir}/snpEff.build.config"
mv "${tmpdir}" "${ROOT}/0_tools/snpeff_Hco_WBPS18"
'
        )"

        echo "${submission}"
        jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
        wait_for_job "${jid}"
    fi

    jobs=()

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue

        mkdir -p "${STAGE}/per_sample/${id}"

        submission="$(
            ROOT="${ROOT}" STAGE="${STAGE}" REFERENCE="${REFERENCE}" ID="${id}" \
            bsub \
                -J "VAR_${id}" \
                -n 4 \
                -R "select[mem>=16000] rusage[mem=16000] span[hosts=1]" \
                -M 16000 \
                -o "${STAGE}/per_sample/${id}/VAR_${id}.%J.out" \
                -e "${STAGE}/per_sample/${id}/VAR_${id}.%J.err" \
                bash -lc '
set -euo pipefail

filtered_tmp="$(mktemp "${STAGE}/filtered_bam/.${ID}.q21.bam.tmp.XXXXXX")"
rm -f "${filtered_tmp}"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/samtools_1.22.sif" \
    samtools view \
    --threads 4 \
    -b \
    -q 21 \
    -o "${filtered_tmp}" \
    "${ROOT}/2_mapping/outputs/${ID}/${ID}.bam"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/samtools_1.22.sif" \
    samtools index \
    -@ 4 \
    -o "${filtered_tmp}.bai" \
    "${filtered_tmp}"

mv "${filtered_tmp}" "${STAGE}/filtered_bam/${ID}.q21.bam"
mv "${filtered_tmp}.bai" "${STAGE}/filtered_bam/${ID}.q21.bam.bai"

mpileup_tmp="$(mktemp "${STAGE}/per_sample/${ID}/.${ID}.mpileup.bcf.tmp.XXXXXX")"
rm -f "${mpileup_tmp}"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/bcftools_1.23.sif" \
    bcftools mpileup \
    --threads 4 \
    -f "${REFERENCE}" \
    -q 21 \
    -a FORMAT/AD,FORMAT/DP,FORMAT/ADF,FORMAT/ADR \
    -Ob \
    -o "${mpileup_tmp}" \
    "${STAGE}/filtered_bam/${ID}.q21.bam"

called_tmp="$(mktemp "${STAGE}/per_sample/${ID}/.${ID}.called.bcf.tmp.XXXXXX")"
rm -f "${called_tmp}"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/bcftools_1.23.sif" \
    bcftools call \
    --threads 4 \
    -m \
    -A \
    -Ob \
    -o "${called_tmp}" \
    "${mpileup_tmp}"

norm_preheader_tmp="$(mktemp "${STAGE}/per_sample/${ID}/.${ID}.norm.preheader.vcf.gz.tmp.XXXXXX")"
rm -f "${norm_preheader_tmp}"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/bcftools_1.23.sif" \
    bcftools norm \
    --threads 4 \
    -f "${REFERENCE}" \
    -m -any \
    -Oz \
    -o "${norm_preheader_tmp}" \
    "${called_tmp}"

rm -f "${mpileup_tmp}" "${called_tmp}"

sample_name_tmp="$(mktemp "${STAGE}/per_sample/${ID}/.${ID}.sample_name.tmp.XXXXXX")"
printf "%s\n" "${ID}" > "${sample_name_tmp}"

norm_tmp="$(mktemp "${STAGE}/per_sample/${ID}/.${ID}.norm.vcf.gz.tmp.XXXXXX")"
rm -f "${norm_tmp}"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/bcftools_1.23.sif" \
    bcftools reheader \
    -s "${sample_name_tmp}" \
    -o "${norm_tmp}" \
    "${norm_preheader_tmp}"

rm -f "${sample_name_tmp}" "${norm_preheader_tmp}"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/bcftools_1.23.sif" \
    bcftools index \
    -f \
    -t \
    "${norm_tmp}"

mv "${norm_tmp}" "${STAGE}/per_sample/${ID}/${ID}.norm.vcf.gz"
mv "${norm_tmp}.tbi" "${STAGE}/per_sample/${ID}/${ID}.norm.vcf.gz.tbi"

annotated_vcf_tmp="$(mktemp "${STAGE}/per_sample/${ID}/.${ID}.annotated.vcf.tmp.XXXXXX")"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/snpeff_5.1d.sif" \
    snpEff ann \
    -noStats \
    -c "${ROOT}/0_tools/snpeff_Hco_WBPS18/snpEff.config" \
    -v Hco_WBPS18 \
    "${STAGE}/per_sample/${ID}/${ID}.norm.vcf.gz" \
    > "${annotated_vcf_tmp}"

annotated_tmp="$(mktemp "${STAGE}/per_sample/${ID}/.${ID}.annotated.vcf.gz.tmp.XXXXXX")"
rm -f "${annotated_tmp}"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/bcftools_1.23.sif" \
    bcftools view \
    -Oz \
    -o "${annotated_tmp}" \
    "${annotated_vcf_tmp}"

rm -f "${annotated_vcf_tmp}"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/bcftools_1.23.sif" \
    bcftools index \
    -f \
    -t \
    "${annotated_tmp}"

mv "${annotated_tmp}" "${STAGE}/per_sample/${ID}/${ID}.annotated.vcf.gz"
mv "${annotated_tmp}.tbi" "${STAGE}/per_sample/${ID}/${ID}.annotated.vcf.gz.tbi"
'
        )"

        echo "${submission}"
        jobs+=("$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")")
    done < "${MANIFEST}"

    for jid in "${jobs[@]}"; do
        wait_for_job "${jid}"
    done

    mkdir -p "${STAGE}/merged"

    submission="$(
        ROOT="${ROOT}" STAGE="${STAGE}" MANIFEST="${MANIFEST}" \
        bsub \
            -J "HcoMergeVCF" \
            -n 4 \
            -R "select[mem>=16000] rusage[mem=16000] span[hosts=1]" \
            -M 16000 \
            -o "${STAGE}/merged/HcoMergeVCF.%J.out" \
            -e "${STAGE}/merged/HcoMergeVCF.%J.err" \
            bash -lc '
set -euo pipefail

vcf_list="$(mktemp "${STAGE}/merged/.annotated_vcfs.list.tmp.XXXXXX")"

while IFS=, read -r id r1 r2; do
    id="${id//$'"'"'\r'"'"'/}"
    [[ "${id}" == "ID" || -z "${id}" ]] && continue
    printf "%s\n" "${STAGE}/per_sample/${id}/${id}.annotated.vcf.gz" >> "${vcf_list}"
done < "${MANIFEST}"

merged_tmp="$(mktemp "${STAGE}/merged/.Hco.multisample.annotated.vcf.gz.tmp.XXXXXX")"
rm -f "${merged_tmp}"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/bcftools_1.23.sif" \
    bcftools merge \
    --threads 4 \
    -m none \
    -l "${vcf_list}" \
    -Oz \
    -o "${merged_tmp}"

rm -f "${vcf_list}"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/bcftools_1.23.sif" \
    bcftools index \
    -f \
    -t \
    "${merged_tmp}"

mv "${merged_tmp}" "${STAGE}/merged/Hco.multisample.annotated.vcf.gz"
mv "${merged_tmp}.tbi" "${STAGE}/merged/Hco.multisample.annotated.vcf.gz.tbi"
'
    )"

    echo "${submission}"
    jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
    wait_for_job "${jid}"

    mv "${STAGE}" "${ROOT}/3_variants"
fi


# ==============================================================================
# 4_FST — pairwise genetic differentiation with Grenedalf v0.6.3 in
# overlapping 5 kb windows with a 2.5 kb step
# ==============================================================================
if [[ ! -d "${ROOT}/4_FST" ]]; then
    echo "Creating and running 4_FST"

    STAGE="$(mktemp -d "${ROOT}/4_FST.tmp.XXXXXX")"

    mkdir -p "${ROOT}/0_tools"

    if [[ ! -f "${ROOT}/0_tools/grenedalf_0.6.3.sif" ]]; then
        grenedalf_tmp="$(find "${ROOT}/0_tools" -maxdepth 1 -type f -name '.grenedalf_0.6.3.sif.tmp.*' -print -quit)"
        if [[ -z "${grenedalf_tmp}" ]]; then
            grenedalf_tmp="$(mktemp "${ROOT}/0_tools/.grenedalf_0.6.3.sif.tmp.XXXXXX")"
        fi
        wget -c \
            -O "${grenedalf_tmp}" \
            "https://depot.galaxyproject.org/singularity/grenedalf%3A0.6.3--hbefcdb2_0"
        mv "${grenedalf_tmp}" "${ROOT}/0_tools/grenedalf_0.6.3.sif"
    fi

    submission="$(
        ROOT="${ROOT}" STAGE="${STAGE}" MANIFEST="${MANIFEST}" \
        bsub \
            -J "HcoFST" \
            -n 16 \
            -R "select[mem>=32000] rusage[mem=32000] span[hosts=1]" \
            -M 32000 \
            -o "${STAGE}/HcoFST.%J.out" \
            -e "${STAGE}/HcoFST.%J.err" \
            bash -lc '
set -euo pipefail

bams=()

while IFS=, read -r id r1 r2; do
    id="${id//$'"'"'\r'"'"'/}"
    [[ "${id}" == "ID" || -z "${id}" ]] && continue
    bams+=("${ROOT}/2_mapping/outputs/${id}/${id}.bam")
done < "${MANIFEST}"

tmpdir="$(mktemp -d "${STAGE}/.windows_5kb.tmp.XXXXXX")"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/grenedalf_0.6.3.sif" \
    grenedalf fst \
    --sam-path "${bams[@]}" \
    --method unbiased-nei \
    --sam-min-map-qual 30 \
    --sam-min-base-qual 30 \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 30 \
    --filter-sample-max-read-depth 300 \
    --pool-sizes 1000 \
    --window-type interval \
    --window-interval-width 5000 \
    --window-interval-stride 2500 \
    --window-average-policy valid-loci \
    --separator-char tab \
    --na-entry nan \
    --threads 16 \
    --compress \
    --out-dir "${tmpdir}"

mv "${tmpdir}" "${STAGE}/windows_5kb"
'
    )"

    echo "${submission}"
    jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
    wait_for_job "${jid}"

    mv "${STAGE}" "${ROOT}/4_FST"
fi

echo "Pipeline completed."
