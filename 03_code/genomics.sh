#!/usr/bin/env bash
#
# PREREQUISITES
#
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
# - if an analysis-stage directory is absent, it is created and launched,
# - if the the run fails, manually remove the last stage directory (the stage failed), and rerun,
# - each output is first written to a .tmp.XXXXXX path and is renamed only after successful completion, so remaining of .tmp.XXXXXX files means that the stage is corrupted.

set -euo pipefail
IFS=$'\n\t'

ROOT="$(pwd -P)"
MANIFEST="$(find "${ROOT}" -maxdepth 1 -type f -name '*.manifest' -print -quit)"
REFERENCE="$(find "${ROOT}/0_refseq" -maxdepth 1 -type f -name '*.fa' -print -quit)"
ANNOTATION="$(find "${ROOT}/0_refseq" -maxdepth 1 -type f -name '*.gtf' -print -quit)"


# ==============================================================================
# 0_samples — sequential acquisition of paired-end FASTQ files with wget -c
# ==============================================================================
if [[ ! -d "${ROOT}/0_samples" ]]; then
    echo "Creating and running 0_samples"
    mkdir -p "${ROOT}/0_samples"

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        r1="${r1//$'\r'/}"
        r2="${r2//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue

        mkdir -p "${ROOT}/0_samples/${id}"

        submission="$(
            ROOT="${ROOT}" ID="${id}" URL_R1="${r1}" URL_R2="${r2}" \
            bsub \
                -J "DL_${id}" \
                -n 1 \
                -R "rusage[mem=2000]" \
                -M 2000 \
                bash -lc '
set -euo pipefail

final_r1="${ROOT}/0_samples/${ID}/${ID}_R1.fastq.gz"
tmp_r1="$(find "${ROOT}/0_samples/${ID}" -maxdepth 1 -type f -name ".${ID}_R1.fastq.gz.tmp.*" -print -quit)"
if [[ -z "${tmp_r1}" ]]; then
    tmp_r1="$(mktemp "${ROOT}/0_samples/${ID}/.${ID}_R1.fastq.gz.tmp.XXXXXX")"
fi
wget -c -O "${tmp_r1}" "${URL_R1}"
mv "${tmp_r1}" "${final_r1}"

final_r2="${ROOT}/0_samples/${ID}/${ID}_R2.fastq.gz"
tmp_r2="$(find "${ROOT}/0_samples/${ID}" -maxdepth 1 -type f -name ".${ID}_R2.fastq.gz.tmp.*" -print -quit)"
if [[ -z "${tmp_r2}" ]]; then
    tmp_r2="$(mktemp "${ROOT}/0_samples/${ID}/.${ID}_R2.fastq.gz.tmp.XXXXXX")"
fi
wget -c -O "${tmp_r2}" "${URL_R2}"
mv "${tmp_r2}" "${final_r2}"
'
        )"

        echo "${submission}"
        jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
        bwait -w "ended(${jid})"
    done < "${MANIFEST}"
fi


# ==============================================================================
# 1_QC — parallel FastQC assessment followed by MultiQC summarisation
# ==============================================================================
if [[ ! -d "${ROOT}/1_QC" ]]; then
    echo "Creating and running 1_QC"

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

    mkdir -p "${ROOT}/1_QC/fastqc"
    jobs=()

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue

        submission="$(
            ROOT="${ROOT}" ID="${id}" \
            bsub \
                -J "QC_${id}" \
                -n 4 \
                -R "rusage[mem=4000]" \
                -M 4000 \
                bash -lc '
set -euo pipefail

tmpdir="$(mktemp -d "${ROOT}/1_QC/fastqc/.${ID}.tmp.XXXXXX")"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/fastqc_0.12.1.sif" \
    fastqc \
    --threads 4 \
    --outdir "${tmpdir}" \
    "${ROOT}/0_samples/${ID}/${ID}_R1.fastq.gz" \
    "${ROOT}/0_samples/${ID}/${ID}_R2.fastq.gz"

mv "${tmpdir}" "${ROOT}/1_QC/fastqc/${ID}"
'
        )"

        echo "${submission}"
        jobs+=("$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")")
    done < "${MANIFEST}"

    for jid in "${jobs[@]}"; do
        bwait -w "ended(${jid})"
    done

    submission="$(
        ROOT="${ROOT}" \
        bsub \
            -J "QC_MultiQC" \
            -n 1 \
            -R "rusage[mem=8000]" \
            -M 8000 \
            bash -lc '
set -euo pipefail

tmpdir="$(mktemp -d "${ROOT}/1_QC/.multiqc.tmp.XXXXXX")"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/multiqc_1.17.sif" \
    multiqc \
    --filename raw_reads_multiqc.html \
    --outdir "${tmpdir}" \
    "${ROOT}/1_QC/fastqc"

mv "${tmpdir}" "${ROOT}/1_QC/multiqc"
'
    )"

    echo "${submission}"
    jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
    bwait -w "ended(${jid})"
fi


# ==============================================================================
# 2_mapping — mapping-helminth v1.0.8 followed by MultiQC summarisation of
# SAMtools flagstat and stats outputs generated by the mapping pipeline
# ==============================================================================
if [[ ! -d "${ROOT}/2_mapping" ]]; then
    echo "Creating and running 2_mapping"

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

    mkdir -p "${ROOT}/2_mapping"

    manifest_tmp="$(mktemp "${ROOT}/2_mapping/.fqpaths.manifest.tmp.XXXXXX")"
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

    mv "${manifest_tmp}" "${ROOT}/2_mapping/fqpaths.manifest"

    submission="$(
        ROOT="${ROOT}" REFERENCE="${REFERENCE}" \
        bsub \
            -J "HcoMap" \
            -n 1 \
            -R "rusage[mem=80000]" \
            -M 80000 \
            bash -lc '
set -euo pipefail

source /etc/profile.d/modules.sh
module load mapping-helminth/v1.0.8

tmpdir="$(mktemp -d "${ROOT}/2_mapping/.outputs.tmp.XXXXXX")"
export NXF_WORK="${ROOT}/2_mapping/work"
export NXF_LOG_FILE=/dev/null
mkdir -p "${NXF_WORK}" "${ROOT}/2_mapping/index_cache"

mapping-helminth \
    --reference "${REFERENCE}" \
    --input "${ROOT}/2_mapping/fqpaths.manifest" \
    --outdir "${tmpdir}" \
    --index_cache "${ROOT}/2_mapping/index_cache" \
    --keep_unmapped true

mv "${tmpdir}" "${ROOT}/2_mapping/outputs"
'
    )"

    echo "${submission}"
    jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
    bwait -w "ended(${jid})"

    submission="$(
        ROOT="${ROOT}" \
        bsub \
            -J "Map_MultiQC" \
            -n 1 \
            -R "rusage[mem=8000]" \
            -M 8000 \
            bash -lc '
set -euo pipefail

tmpdir="$(mktemp -d "${ROOT}/2_mapping/.multiqc.tmp.XXXXXX")"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/multiqc_1.17.sif" \
    multiqc \
    --filename mapping_multiqc.html \
    --outdir "${tmpdir}" \
    "${ROOT}/2_mapping/outputs"

mv "${tmpdir}" "${ROOT}/2_mapping/multiqc"
'
    )"

    echo "${submission}"
    jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
    bwait -w "ended(${jid})"
fi


# ==============================================================================
# 3_variants — strict MAPQ >20 filtering, BCFtools v1.23 pile-up generation,
# multiallelic calling, biallelic normalisation and SnpEff v5.1d annotation
# SNPs and short insertions/deletions are retained.
# ==============================================================================
if [[ ! -d "${ROOT}/3_variants" ]]; then
    echo "Creating and running 3_variants"

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

    mkdir -p "${ROOT}/3_variants/filtered_bam"
    mkdir -p "${ROOT}/3_variants/per_sample"

    # The species-specific database is built from the same WBPS18 genome and
    # annotation used for mapping and subsequent genomic interpretation.
    if [[ ! -d "${ROOT}/0_tools/snpeff_Hco_WBPS18" ]]; then
        submission="$(
            ROOT="${ROOT}" REFERENCE="${REFERENCE}" ANNOTATION="${ANNOTATION}" \
            bsub \
                -J "HcoSnpEffDB" \
                -n 4 \
                -R "rusage[mem=12000]" \
                -M 12000 \
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
        bwait -w "ended(${jid})"
    fi

    jobs=()

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue

        mkdir -p "${ROOT}/3_variants/per_sample/${id}"

        submission="$(
            ROOT="${ROOT}" REFERENCE="${REFERENCE}" ID="${id}" \
            bsub \
                -J "VAR_${id}" \
                -n 4 \
                -R "rusage[mem=16000]" \
                -M 16000 \
                bash -lc '
set -euo pipefail

filtered_tmp="$(mktemp "${ROOT}/3_variants/filtered_bam/.${ID}.q21.bam.tmp.XXXXXX")"
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

mv "${filtered_tmp}" "${ROOT}/3_variants/filtered_bam/${ID}.q21.bam"
mv "${filtered_tmp}.bai" "${ROOT}/3_variants/filtered_bam/${ID}.q21.bam.bai"

mpileup_tmp="$(mktemp "${ROOT}/3_variants/per_sample/${ID}/.${ID}.mpileup.bcf.tmp.XXXXXX")"
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
    "${ROOT}/3_variants/filtered_bam/${ID}.q21.bam"

called_tmp="$(mktemp "${ROOT}/3_variants/per_sample/${ID}/.${ID}.called.bcf.tmp.XXXXXX")"
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

norm_preheader_tmp="$(mktemp "${ROOT}/3_variants/per_sample/${ID}/.${ID}.norm.preheader.vcf.gz.tmp.XXXXXX")"
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

sample_name_tmp="$(mktemp "${ROOT}/3_variants/per_sample/${ID}/.${ID}.sample_name.tmp.XXXXXX")"
printf "%s\n" "${ID}" > "${sample_name_tmp}"

norm_tmp="$(mktemp "${ROOT}/3_variants/per_sample/${ID}/.${ID}.norm.vcf.gz.tmp.XXXXXX")"
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

mv "${norm_tmp}" "${ROOT}/3_variants/per_sample/${ID}/${ID}.norm.vcf.gz"
mv "${norm_tmp}.tbi" "${ROOT}/3_variants/per_sample/${ID}/${ID}.norm.vcf.gz.tbi"

annotated_vcf_tmp="$(mktemp "${ROOT}/3_variants/per_sample/${ID}/.${ID}.annotated.vcf.tmp.XXXXXX")"

singularity exec \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/snpeff_5.1d.sif" \
    snpEff ann \
    -noStats \
    -c "${ROOT}/0_tools/snpeff_Hco_WBPS18/snpEff.config" \
    -v Hco_WBPS18 \
    "${ROOT}/3_variants/per_sample/${ID}/${ID}.norm.vcf.gz" \
    > "${annotated_vcf_tmp}"

annotated_tmp="$(mktemp "${ROOT}/3_variants/per_sample/${ID}/.${ID}.annotated.vcf.gz.tmp.XXXXXX")"
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

mv "${annotated_tmp}" "${ROOT}/3_variants/per_sample/${ID}/${ID}.annotated.vcf.gz"
mv "${annotated_tmp}.tbi" "${ROOT}/3_variants/per_sample/${ID}/${ID}.annotated.vcf.gz.tbi"
'
        )"

        echo "${submission}"
        jobs+=("$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")")
    done < "${MANIFEST}"

    for jid in "${jobs[@]}"; do
        bwait -w "ended(${jid})"
    done

    mkdir -p "${ROOT}/3_variants/merged"

    submission="$(
        ROOT="${ROOT}" MANIFEST="${MANIFEST}" \
        bsub \
            -J "HcoMergeVCF" \
            -n 4 \
            -R "rusage[mem=16000]" \
            -M 16000 \
            bash -lc '
set -euo pipefail

vcf_list="$(mktemp "${ROOT}/3_variants/merged/.annotated_vcfs.list.tmp.XXXXXX")"

while IFS=, read -r id r1 r2; do
    id="${id//$'"'"'\r'"'"'/}"
    [[ "${id}" == "ID" || -z "${id}" ]] && continue
    printf "%s\n" "${ROOT}/3_variants/per_sample/${id}/${id}.annotated.vcf.gz" >> "${vcf_list}"
done < "${MANIFEST}"

merged_tmp="$(mktemp "${ROOT}/3_variants/merged/.Hco.multisample.annotated.vcf.gz.tmp.XXXXXX")"
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

mv "${merged_tmp}" "${ROOT}/3_variants/merged/Hco.multisample.annotated.vcf.gz"
mv "${merged_tmp}.tbi" "${ROOT}/3_variants/merged/Hco.multisample.annotated.vcf.gz.tbi"
'
    )"

    echo "${submission}"
    jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
    bwait -w "ended(${jid})"
fi


# ==============================================================================
# 4_FST — pairwise genetic differentiation with Grenedalf v0.6.3 in
# overlapping 5 kb windows with a 2.5 kb step
# ==============================================================================
if [[ ! -d "${ROOT}/4_FST" ]]; then
    echo "Creating and running 4_FST"

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

    mkdir -p "${ROOT}/4_FST"

    submission="$(
        ROOT="${ROOT}" MANIFEST="${MANIFEST}" \
        bsub \
            -J "HcoFST" \
            -n 16 \
            -R "rusage[mem=32000]" \
            -M 32000 \
            bash -lc '
set -euo pipefail

bams=()

while IFS=, read -r id r1 r2; do
    id="${id//$'"'"'\r'"'"'/}"
    [[ "${id}" == "ID" || -z "${id}" ]] && continue
    bams+=("${ROOT}/2_mapping/outputs/${id}/${id}.bam")
done < "${MANIFEST}"

tmpdir="$(mktemp -d "${ROOT}/4_FST/.windows_5kb.tmp.XXXXXX")"

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

mv "${tmpdir}" "${ROOT}/4_FST/windows_5kb"
'
    )"

    echo "${submission}"
    jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
    bwait -w "ended(${jid})"
fi

echo "Pipeline completed."
