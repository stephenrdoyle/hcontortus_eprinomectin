#!/usr/bin/env bash
#
# TRANSCRIPTOMICS PIPELINE — Haemonchus contortus RNA-seq
#
# This script consolidates the supplied draft pipeline (main.sh + 01_fastqc.sh,
# 02_multiqc.sh, 03_cutadapt.sh, 04_hisat2.sh, 05_MarkDuplicates.sh,
# 05_quality.sh and 11_featureCounts.sh/R) into one LSF orchestrator.
#
# SCOPE
# - raw paired-end FASTQ acquisition from transcriptomics.csv.manifest
# - FastQC 0.12.1
# - MultiQC 1.14 (version retained from the supplied draft)
# - Cutadapt 4.3 with the supplied Illumina adapter sequences and otherwise
#   default trimming settings
# - HISAT2 2.2.1 mapping to H. contortus PRJEB506 / WBPS18
# - SAM -> BAM conversion, coordinate sorting and indexing with SAMtools 1.19
# - duplicate marking with Sambamba 1.0.1
# - mapping QC with SAMtools 1.19 flagstat and stats
# - paired-end gene counting with featureCounts from Rsubread 2.16.1
#
# IMPORTANT METHOD BOUNDARY
# The supplied draft ends at featureCounts. The later article Methods sections
# (DESeq2 normalisation, sample concordance, BUN-M3 exclusion, DGE and enrichment)
# are therefore not invented here without their original analysis scripts.
# In particular, Hc_T_R_BUN_M_3 is retained through the raw-count table and is
# only to be excluded at the downstream differential-expression stage, exactly
# as described in the Methods.
#
# PREREQUISITES
# - HPC cluster using LSF: bsub, bwait and bjobs available
# - Singularity or Apptainer available in the launch environment
# - wget, gzip and standard GNU/POSIX shell utilities available
# - current launch directory contains:
#     ./transcriptomics.csv.manifest
#     ./0_refseq/<one WBPS18 PRJEB506 genome .fa>
#     ./0_refseq/<one WBPS18 PRJEB506 canonical geneset .gtf>
#
# MANIFEST FORMAT
# ID,R1,R2
# Hc_T_R_ARA_F_1,ftp.sra.ebi.ac.uk/..._1.fastq.gz,ftp.sra.ebi.ac.uk/..._2.fastq.gz
#
# URLs without an explicit scheme are interpreted as ftp:// URLs. This matches
# the supplied transcriptomics.csv.manifest.
#
# EXECUTION PRINCIPLE
# - each stage is created as <stage>.tmp.XXXXXX;
# - the stage is renamed to its final directory only after every LSF job in that
#   stage has completed successfully and expected outputs have been validated;
# - a pre-existing final stage directory is treated as completed and skipped;
# - a stale temporary stage directory causes an explicit stop: remove it after
#   inspection before rerunning;
# - sample identity always comes from the manifest ID field. Filenames are never
#   truncated at underscores or reconstructed heuristically;
# - all stage transitions use deterministic paths and validate required files.
#
# RECOMMENDED LAUNCH
# bsub -q long -J transcriptomics -n 1 \
#   -R "select[mem>=2000] rusage[mem=2000] span[hosts=1]" -M 2000 \
#   -o "transcriptomics.out" -e "transcriptomics.err" bash transcriptomics.sh

set -euo pipefail
IFS=$'\n\t'

ROOT="$(pwd -P)"
MANIFEST="${ROOT}/transcriptomics.csv.manifest"
REFDIR="${ROOT}/0_refseq"
TOOLDIR="${ROOT}/0_tools"

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_nonempty_file() {
    [[ -s "$1" ]] || die "missing or empty file: $1"
}

begin_stage() {
    local name="$1"
    local final="${ROOT}/${name}"
    local stale

    if [[ -d "${final}" ]]; then
        echo "Skipping ${name}: final stage directory already exists."
        return 1
    fi

    stale="$(find "${ROOT}" -maxdepth 1 -type d -name "${name}.tmp.*" -print -quit)"
    if [[ -n "${stale}" ]]; then
        die "stale temporary stage detected: ${stale}. Inspect/remove it before rerunning."
    fi

    return 0
}

wait_for_job() {
    local jid="$1"
    local status

    [[ "${jid}" =~ ^[0-9]+$ ]] || die "invalid or missing LSF job ID: '${jid}'"

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
                die "LSF job ${jid} ended with status EXIT. Pipeline stopped; inspect its .err/.out files and the temporary stage."
                ;;
            *)
                sleep 2
                ;;
        esac
    done
}

extract_job_id() {
    local submission="$1"
    local jid
    jid="$(sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' <<< "${submission}")"
    [[ "${jid}" =~ ^[0-9]+$ ]] || die "could not parse LSF job ID from: ${submission}"
    printf '%s\n' "${jid}"
}

ensure_sif() {
    local filename="$1"
    local galaxy_url="$2"
    local docker_uri="$3"
    local final="${TOOLDIR}/${filename}"
    local tmp

    mkdir -p "${TOOLDIR}"

    if [[ -s "${final}" ]]; then
        return 0
    fi

    tmp="$(find "${TOOLDIR}" -maxdepth 1 -type f -name ".${filename}.tmp.*" -print -quit)"
    if [[ -z "${tmp}" ]]; then
        tmp="$(mktemp "${TOOLDIR}/.${filename}.tmp.XXXXXX")"
    fi

    echo "Obtaining ${filename}"

    if ! wget -c -O "${tmp}" "${galaxy_url}"; then
        echo "Galaxy Singularity download failed; falling back to ${docker_uri}" >&2
        rm -f "${tmp}"
        tmp="$(mktemp "${TOOLDIR}/.${filename}.tmp.XXXXXX")"
        rm -f "${tmp}"
        "${CONTAINER_RUNTIME}" pull "${tmp}" "${docker_uri}"
    fi

    [[ -s "${tmp}" ]] || die "container download produced an empty file: ${tmp}"
    mv "${tmp}" "${final}"
}

# ------------------------------------------------------------------------------
# Global validation
# ------------------------------------------------------------------------------
command -v bsub >/dev/null 2>&1 || die "bsub is not available"
command -v bwait >/dev/null 2>&1 || die "bwait is not available"
command -v bjobs >/dev/null 2>&1 || die "bjobs is not available"
command -v wget >/dev/null 2>&1 || die "wget is not available"
command -v gzip >/dev/null 2>&1 || die "gzip is not available"

if command -v singularity >/dev/null 2>&1; then
    CONTAINER_RUNTIME="$(command -v singularity)"
elif command -v apptainer >/dev/null 2>&1; then
    CONTAINER_RUNTIME="$(command -v apptainer)"
else
    die "neither singularity nor apptainer is available; load the relevant module before launching"
fi

require_nonempty_file "${MANIFEST}"
[[ -d "${REFDIR}" ]] || die "reference directory not found: ${REFDIR}"

header="$(head -n 1 "${MANIFEST}" | tr -d '\r')"
[[ "${header}" == "ID,R1,R2" ]] || die "manifest header must be exactly: ID,R1,R2"

manifest_problem="$(
    awk -F',' '
        NR == 1 { next }
        {
            gsub(/\r/, "", $0)
            if ($0 == "") next
            if (NF != 3 || $1 == "" || $2 == "" || $3 == "") {
                print "malformed line " NR ": " $0
                exit
            }
            if ($1 !~ /^[A-Za-z0-9._-]+$/) {
                print "unsafe sample ID on line " NR ": " $1
                exit
            }
            if (seen[$1]++) {
                print "duplicate sample ID on line " NR ": " $1
                exit
            }
        }
    ' "${MANIFEST}"
)"
[[ -z "${manifest_problem}" ]] || die "manifest validation failed: ${manifest_problem}"

SAMPLE_COUNT="$(awk -F',' 'NR > 1 {gsub(/\r/, "", $0); if ($1 != "") n++} END {print n+0}' "${MANIFEST}")"
(( SAMPLE_COUNT > 0 )) || die "manifest contains no samples"

mapfile -t REFERENCES < <(find "${REFDIR}" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fasta' \) -print | sort)
mapfile -t ANNOTATIONS < <(find "${REFDIR}" -maxdepth 1 -type f -name '*.gtf' -print | sort)

(( ${#REFERENCES[@]} == 1 )) || die "expected exactly one reference .fa/.fasta in ${REFDIR}; found ${#REFERENCES[@]}"
(( ${#ANNOTATIONS[@]} == 1 )) || die "expected exactly one annotation .gtf in ${REFDIR}; found ${#ANNOTATIONS[@]}"

REFERENCE="${REFERENCES[0]}"
ANNOTATION="${ANNOTATIONS[0]}"

ref_base="$(basename "${REFERENCE}")"
gtf_base="$(basename "${ANNOTATION}")"
[[ "${ref_base}" == *PRJEB506* && "${ref_base}" == *WBPS18* ]] \
    || die "reference must be H. contortus PRJEB506 WBPS18; found ${ref_base}"
[[ "${gtf_base}" == *PRJEB506* && "${gtf_base}" == *WBPS18* ]] \
    || die "annotation must be H. contortus PRJEB506 WBPS18; found ${gtf_base}"

mkdir -p "${TOOLDIR}"

# Exact, immutable BioContainers used by the supplied Methods/draft versions.
FASTQC_SIF="${TOOLDIR}/fastqc_0.12.1.sif"
MULTIQC_SIF="${TOOLDIR}/multiqc_1.14.sif"
CUTADAPT_SIF="${TOOLDIR}/cutadapt_4.3.sif"
HISAT2_SIF="${TOOLDIR}/hisat2_2.2.1.sif"
SAMTOOLS_SIF="${TOOLDIR}/samtools_1.19.sif"
SAMBAMBA_SIF="${TOOLDIR}/sambamba_1.0.1.sif"
RSUBREAD_SIF="${TOOLDIR}/rsubread_2.16.1.sif"

# ==============================================================================
# 0_samples — sequential acquisition from transcriptomics.csv.manifest
# ==============================================================================
if begin_stage "0_samples"; then
    echo "Creating and running 0_samples (${SAMPLE_COUNT} samples)"
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

normalise_url() {
    case "$1" in
        *://*) printf "%s\n" "$1" ;;
        *)     printf "ftp://%s\n" "$1" ;;
    esac
}

u1="$(normalise_url "${URL_R1}")"
u2="$(normalise_url "${URL_R2}")"

final_r1="${STAGE}/${ID}/${ID}_R1.fastq.gz"
tmp_r1="$(find "${STAGE}/${ID}" -maxdepth 1 -type f -name ".${ID}_R1.fastq.gz.tmp.*" -print -quit)"
if [[ -z "${tmp_r1}" ]]; then
    tmp_r1="$(mktemp "${STAGE}/${ID}/.${ID}_R1.fastq.gz.tmp.XXXXXX")"
fi
wget -c -O "${tmp_r1}" "${u1}"
gzip -t "${tmp_r1}"
[[ -s "${tmp_r1}" ]]
mv "${tmp_r1}" "${final_r1}"

final_r2="${STAGE}/${ID}/${ID}_R2.fastq.gz"
tmp_r2="$(find "${STAGE}/${ID}" -maxdepth 1 -type f -name ".${ID}_R2.fastq.gz.tmp.*" -print -quit)"
if [[ -z "${tmp_r2}" ]]; then
    tmp_r2="$(mktemp "${STAGE}/${ID}/.${ID}_R2.fastq.gz.tmp.XXXXXX")"
fi
wget -c -O "${tmp_r2}" "${u2}"
gzip -t "${tmp_r2}"
[[ -s "${tmp_r2}" ]]
mv "${tmp_r2}" "${final_r2}"
'
        )"

        echo "${submission}"
        jid="$(extract_job_id "${submission}")"
        wait_for_job "${jid}"
    done < "${MANIFEST}"

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue
        require_nonempty_file "${STAGE}/${id}/${id}_R1.fastq.gz"
        require_nonempty_file "${STAGE}/${id}/${id}_R2.fastq.gz"
        gzip -t "${STAGE}/${id}/${id}_R1.fastq.gz"
        gzip -t "${STAGE}/${id}/${id}_R2.fastq.gz"
    done < "${MANIFEST}"

    mv "${STAGE}" "${ROOT}/0_samples"
fi

# ==============================================================================
# 01_fastqc — raw-read QC with FastQC 0.12.1
# ==============================================================================
if begin_stage "01_fastqc"; then
    echo "Creating and running 01_fastqc"
    STAGE="$(mktemp -d "${ROOT}/01_fastqc.tmp.XXXXXX")"

    ensure_sif \
        "fastqc_0.12.1.sif" \
        "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0" \
        "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    jobs=()

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue
        mkdir -p "${STAGE}/${id}"

        submission="$(
            ROOT="${ROOT}" STAGE="${STAGE}" ID="${id}" RUNTIME="${CONTAINER_RUNTIME}" FASTQC_SIF="${FASTQC_SIF}" \
            bsub \
                -J "FQC_${id}" \
                -n 4 \
                -R "select[mem>=4000] rusage[mem=4000] span[hosts=1]" \
                -M 4000 \
                -o "${STAGE}/${id}/FQC_${id}.%J.out" \
                -e "${STAGE}/${id}/FQC_${id}.%J.err" \
                bash -lc '
set -euo pipefail

tmpdir="$(mktemp -d "${STAGE}/${ID}/.fastqc.tmp.XXXXXX")"

"${RUNTIME}" exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${FASTQC_SIF}" \
    fastqc \
    --threads 4 \
    --outdir "${tmpdir}" \
    "${ROOT}/0_samples/${ID}/${ID}_R1.fastq.gz" \
    "${ROOT}/0_samples/${ID}/${ID}_R2.fastq.gz"

for mate in R1 R2; do
    [[ -s "${tmpdir}/${ID}_${mate}_fastqc.html" ]]
    [[ -s "${tmpdir}/${ID}_${mate}_fastqc.zip" ]]
done

mv "${tmpdir}" "${STAGE}/${ID}/results"
'
        )"

        echo "${submission}"
        jobs+=("$(extract_job_id "${submission}")")
    done < "${MANIFEST}"

    for jid in "${jobs[@]}"; do
        wait_for_job "${jid}"
    done

    mv "${STAGE}" "${ROOT}/01_fastqc"
fi

# ==============================================================================
# 02_multiqc — aggregate raw FastQC reports with MultiQC 1.14
# ==============================================================================
if begin_stage "02_multiqc"; then
    echo "Creating and running 02_multiqc"
    STAGE="$(mktemp -d "${ROOT}/02_multiqc.tmp.XXXXXX")"

    ensure_sif \
        "multiqc_1.14.sif" \
        "https://depot.galaxyproject.org/singularity/multiqc%3A1.14--pyhdfd78af_0" \
        "docker://quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"

    submission="$(
        ROOT="${ROOT}" STAGE="${STAGE}" RUNTIME="${CONTAINER_RUNTIME}" MULTIQC_SIF="${MULTIQC_SIF}" \
        bsub \
            -J "RNA_MultiQC" \
            -n 1 \
            -R "select[mem>=4000] rusage[mem=4000] span[hosts=1]" \
            -M 4000 \
            -o "${STAGE}/RNA_MultiQC.%J.out" \
            -e "${STAGE}/RNA_MultiQC.%J.err" \
            bash -lc '
set -euo pipefail

tmpdir="$(mktemp -d "${STAGE}/.multiqc.tmp.XXXXXX")"

"${RUNTIME}" exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${MULTIQC_SIF}" \
    multiqc \
    --filename raw_reads_multiqc.html \
    --outdir "${tmpdir}" \
    "${ROOT}/01_fastqc"

[[ -s "${tmpdir}/raw_reads_multiqc.html" ]]
mv "${tmpdir}" "${STAGE}/results"
'
    )"

    echo "${submission}"
    jid="$(extract_job_id "${submission}")"
    wait_for_job "${jid}"

    mv "${STAGE}" "${ROOT}/02_multiqc"
fi

# ==============================================================================
# 03_cutadapt — paired-end adapter trimming with Cutadapt 4.3
# Supplied adapter sequences are retained; all other trimming behaviour is left
# at Cutadapt defaults, matching the Methods wording and the draft script.
# ==============================================================================
if begin_stage "03_cutadapt"; then
    echo "Creating and running 03_cutadapt"
    STAGE="$(mktemp -d "${ROOT}/03_cutadapt.tmp.XXXXXX")"

    ensure_sif \
        "cutadapt_4.3.sif" \
        "https://depot.galaxyproject.org/singularity/cutadapt%3A4.3--py310h1425a21_0" \
        "docker://quay.io/biocontainers/cutadapt:4.3--py310h1425a21_0"

    ADAPTER1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    ADAPTER2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    jobs=()

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue
        mkdir -p "${STAGE}/${id}"

        submission="$(
            ROOT="${ROOT}" STAGE="${STAGE}" ID="${id}" RUNTIME="${CONTAINER_RUNTIME}" CUTADAPT_SIF="${CUTADAPT_SIF}" ADAPTER1="${ADAPTER1}" ADAPTER2="${ADAPTER2}" \
            bsub \
                -J "CUT_${id}" \
                -n 4 \
                -R "select[mem>=4000] rusage[mem=4000] span[hosts=1]" \
                -M 4000 \
                -o "${STAGE}/${id}/CUT_${id}.%J.out" \
                -e "${STAGE}/${id}/CUT_${id}.%J.err" \
                bash -lc '
set -euo pipefail

tmpdir="$(mktemp -d "${STAGE}/${ID}/.cutadapt.tmp.XXXXXX")"
out1="${tmpdir}/${ID}_R1.trim.fastq.gz"
out2="${tmpdir}/${ID}_R2.trim.fastq.gz"

"${RUNTIME}" exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${CUTADAPT_SIF}" \
    cutadapt \
    -j 4 \
    -a "${ADAPTER1}" \
    -A "${ADAPTER2}" \
    -o "${out1}" \
    -p "${out2}" \
    "${ROOT}/0_samples/${ID}/${ID}_R1.fastq.gz" \
    "${ROOT}/0_samples/${ID}/${ID}_R2.fastq.gz"

gzip -t "${out1}"
gzip -t "${out2}"
[[ -s "${out1}" && -s "${out2}" ]]

mv "${out1}" "${STAGE}/${ID}/${ID}_R1.trim.fastq.gz"
mv "${out2}" "${STAGE}/${ID}/${ID}_R2.trim.fastq.gz"
rmdir "${tmpdir}"
'
        )"

        echo "${submission}"
        jobs+=("$(extract_job_id "${submission}")")
    done < "${MANIFEST}"

    for jid in "${jobs[@]}"; do
        wait_for_job "${jid}"
    done

    mv "${STAGE}" "${ROOT}/03_cutadapt"
fi

# ==============================================================================
# 04_hisat2 — HISAT2 2.2.1 index + mapping; SAM -> sorted/indexed BAM with
# SAMtools 1.19. The draft read-group intent is retained using valid HISAT2
# syntax: --rg-id is required for the @RG line, with SM/LB/PL supplied via --rg.
# ==============================================================================
if begin_stage "04_hisat2"; then
    echo "Creating and running 04_hisat2"
    STAGE="$(mktemp -d "${ROOT}/04_hisat2.tmp.XXXXXX")"
    mkdir -p "${STAGE}/index" "${STAGE}/bam" "${STAGE}/logs"

    ensure_sif \
        "hisat2_2.2.1.sif" \
        "https://depot.galaxyproject.org/singularity/hisat2%3A2.2.1--h87f3376_4" \
        "docker://quay.io/biocontainers/hisat2:2.2.1--h87f3376_4"

    ensure_sif \
        "samtools_1.19.sif" \
        "https://depot.galaxyproject.org/singularity/samtools%3A1.19--h50ea8bc_0" \
        "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"

    submission="$(
        ROOT="${ROOT}" STAGE="${STAGE}" REFERENCE="${REFERENCE}" RUNTIME="${CONTAINER_RUNTIME}" HISAT2_SIF="${HISAT2_SIF}" \
        bsub \
            -J "HISAT2_index" \
            -n 8 \
            -R "select[mem>=16000] rusage[mem=16000] span[hosts=1]" \
            -M 16000 \
            -o "${STAGE}/HISAT2_index.%J.out" \
            -e "${STAGE}/HISAT2_index.%J.err" \
            bash -lc '
set -euo pipefail

tmpdir="$(mktemp -d "${STAGE}/index/.build.tmp.XXXXXX")"
prefix="${tmpdir}/hisat2idx"

"${RUNTIME}" exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${HISAT2_SIF}" \
    hisat2-build \
    -p 8 \
    "${REFERENCE}" \
    "${prefix}"

count="$(find "${tmpdir}" -maxdepth 1 -type f -name '"'"'hisat2idx.*.ht2*'"'"' | wc -l)"
(( count >= 8 ))

mv "${tmpdir}" "${STAGE}/index/WBPS18"
'
    )"

    echo "${submission}"
    jid="$(extract_job_id "${submission}")"
    wait_for_job "${jid}"

    jobs=()

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue
        mkdir -p "${STAGE}/logs/${id}"

        submission="$(
            ROOT="${ROOT}" STAGE="${STAGE}" ID="${id}" RUNTIME="${CONTAINER_RUNTIME}" HISAT2_SIF="${HISAT2_SIF}" SAMTOOLS_SIF="${SAMTOOLS_SIF}" \
            bsub \
                -J "MAP_${id}" \
                -n 8 \
                -R "select[mem>=16000] rusage[mem=16000] span[hosts=1]" \
                -M 16000 \
                -o "${STAGE}/logs/${id}/MAP_${id}.%J.out" \
                -e "${STAGE}/logs/${id}/MAP_${id}.%J.err" \
                bash -lc '
set -euo pipefail

tmpdir="$(mktemp -d "${STAGE}/bam/.${ID}.tmp.XXXXXX")"
sam="${tmpdir}/${ID}.sam"
unsorted="${tmpdir}/${ID}.bam"
sorted="${tmpdir}/${ID}.sorted.bam"

"${RUNTIME}" exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${HISAT2_SIF}" \
    hisat2 \
    -p 8 \
    --rg-id "${ID}" \
    --rg "SM:${ID}" \
    --rg "LB:library" \
    --rg "PL:illumina" \
    -x "${STAGE}/index/WBPS18/hisat2idx" \
    -1 "${ROOT}/03_cutadapt/${ID}/${ID}_R1.trim.fastq.gz" \
    -2 "${ROOT}/03_cutadapt/${ID}/${ID}_R2.trim.fastq.gz" \
    -S "${sam}"

[[ -s "${sam}" ]]

"${RUNTIME}" exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${SAMTOOLS_SIF}" \
    samtools view \
    -@ 4 \
    -bS \
    -o "${unsorted}" \
    "${sam}"

[[ -s "${unsorted}" ]]
rm "${sam}"

"${RUNTIME}" exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${SAMTOOLS_SIF}" \
    samtools sort \
    -@ 8 \
    -o "${sorted}" \
    "${unsorted}"

[[ -s "${sorted}" ]]
rm "${unsorted}"

"${RUNTIME}" exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${SAMTOOLS_SIF}" \
    samtools index \
    -@ 8 \
    "${sorted}"

[[ -s "${sorted}.bai" ]]

mv "${sorted}" "${STAGE}/bam/${ID}.sorted.bam"
mv "${sorted}.bai" "${STAGE}/bam/${ID}.sorted.bam.bai"
rmdir "${tmpdir}"
'
        )"

        echo "${submission}"
        jobs+=("$(extract_job_id "${submission}")")
    done < "${MANIFEST}"

    for jid in "${jobs[@]}"; do
        wait_for_job "${jid}"
    done

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue
        require_nonempty_file "${STAGE}/bam/${id}.sorted.bam"
        require_nonempty_file "${STAGE}/bam/${id}.sorted.bam.bai"
    done < "${MANIFEST}"

    mv "${STAGE}" "${ROOT}/04_hisat2"
fi

# ==============================================================================
# 05_MarkDuplicates — duplicate marking with Sambamba 1.0.1
# ==============================================================================
if begin_stage "05_MarkDuplicates"; then
    echo "Creating and running 05_MarkDuplicates"
    STAGE="$(mktemp -d "${ROOT}/05_MarkDuplicates.tmp.XXXXXX")"
    mkdir -p "${STAGE}/bam" "${STAGE}/logs"

    ensure_sif \
        "sambamba_1.0.1.sif" \
        "https://depot.galaxyproject.org/singularity/sambamba%3A1.0.1--h6f6fda4_0" \
        "docker://quay.io/biocontainers/sambamba:1.0.1--h6f6fda4_0"

    ensure_sif \
        "samtools_1.19.sif" \
        "https://depot.galaxyproject.org/singularity/samtools%3A1.19--h50ea8bc_0" \
        "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"

    jobs=()

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue
        mkdir -p "${STAGE}/logs/${id}"

        submission="$(
            ROOT="${ROOT}" STAGE="${STAGE}" ID="${id}" RUNTIME="${CONTAINER_RUNTIME}" SAMBAMBA_SIF="${SAMBAMBA_SIF}" SAMTOOLS_SIF="${SAMTOOLS_SIF}" \
            bsub \
                -J "MKD_${id}" \
                -n 8 \
                -R "select[mem>=16000] rusage[mem=16000] span[hosts=1]" \
                -M 16000 \
                -o "${STAGE}/logs/${id}/MKD_${id}.%J.out" \
                -e "${STAGE}/logs/${id}/MKD_${id}.%J.err" \
                bash -lc '
set -euo pipefail

tmpdir="$(mktemp -d "${STAGE}/bam/.${ID}.tmp.XXXXXX")"
out="${tmpdir}/${ID}.mkdup.sorted.bam"

"${RUNTIME}" exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${SAMBAMBA_SIF}" \
    sambamba markdup \
    -t=8 \
    "${ROOT}/04_hisat2/bam/${ID}.sorted.bam" \
    "${out}"

[[ -s "${out}" ]]

"${RUNTIME}" exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${SAMTOOLS_SIF}" \
    samtools index \
    -@ 8 \
    "${out}"

[[ -s "${out}.bai" ]]

mv "${out}" "${STAGE}/bam/${ID}.mkdup.sorted.bam"
mv "${out}.bai" "${STAGE}/bam/${ID}.mkdup.sorted.bam.bai"
rmdir "${tmpdir}"
'
        )"

        echo "${submission}"
        jobs+=("$(extract_job_id "${submission}")")
    done < "${MANIFEST}"

    for jid in "${jobs[@]}"; do
        wait_for_job "${jid}"
    done

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue
        require_nonempty_file "${STAGE}/bam/${id}.mkdup.sorted.bam"
        require_nonempty_file "${STAGE}/bam/${id}.mkdup.sorted.bam.bai"
    done < "${MANIFEST}"

    mv "${STAGE}" "${ROOT}/05_MarkDuplicates"
fi

# ==============================================================================
# 05_quality — mapping QC with SAMtools 1.19 flagstat and stats
# The raw SAMtools reports are the canonical QC outputs required by the Methods.
# A conservative CSV summary is also produced; it avoids inferring a
# "uniquely-mapped" count from flagstat secondary-alignments, which the draft
# did not actually measure.
# ==============================================================================
if begin_stage "05_quality"; then
    echo "Creating and running 05_quality"
    STAGE="$(mktemp -d "${ROOT}/05_quality.tmp.XXXXXX")"
    mkdir -p "${STAGE}/per_sample"

    ensure_sif \
        "samtools_1.19.sif" \
        "https://depot.galaxyproject.org/singularity/samtools%3A1.19--h50ea8bc_0" \
        "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"

    jobs=()

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue
        mkdir -p "${STAGE}/per_sample/${id}"

        submission="$(
            ROOT="${ROOT}" STAGE="${STAGE}" ID="${id}" RUNTIME="${CONTAINER_RUNTIME}" SAMTOOLS_SIF="${SAMTOOLS_SIF}" \
            bsub \
                -J "QCmap_${id}" \
                -n 4 \
                -R "select[mem>=4000] rusage[mem=4000] span[hosts=1]" \
                -M 4000 \
                -o "${STAGE}/per_sample/${id}/QCmap_${id}.%J.out" \
                -e "${STAGE}/per_sample/${id}/QCmap_${id}.%J.err" \
                bash -lc '
set -euo pipefail

bam="${ROOT}/05_MarkDuplicates/bam/${ID}.mkdup.sorted.bam"
flag_tmp="$(mktemp "${STAGE}/per_sample/${ID}/.${ID}.flagstat.tmp.XXXXXX")"
stats_tmp="$(mktemp "${STAGE}/per_sample/${ID}/.${ID}.stats.tmp.XXXXXX")"

"${RUNTIME}" exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${SAMTOOLS_SIF}" \
    samtools flagstat \
    -@ 4 \
    "${bam}" \
    > "${flag_tmp}"

"${RUNTIME}" exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${SAMTOOLS_SIF}" \
    samtools stats \
    -@ 4 \
    "${bam}" \
    > "${stats_tmp}"

[[ -s "${flag_tmp}" && -s "${stats_tmp}" ]]
mv "${flag_tmp}" "${STAGE}/per_sample/${ID}/${ID}_flagstat.txt"
mv "${stats_tmp}" "${STAGE}/per_sample/${ID}/${ID}_stats.txt"
'
        )"

        echo "${submission}"
        jobs+=("$(extract_job_id "${submission}")")
    done < "${MANIFEST}"

    for jid in "${jobs[@]}"; do
        wait_for_job "${jid}"
    done

    summary_tmp="$(mktemp "${STAGE}/.summary_metrics.csv.tmp.XXXXXX")"
    printf 'Sample,TotalReads,MappedReads,AlignmentRatePct,ProperlyPairedReads,ProperlyPairedPct,DuplicateReads,DuplicationRatePct,BasesMappedCigar\n' > "${summary_tmp}"

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue

        flag="${STAGE}/per_sample/${id}/${id}_flagstat.txt"
        stats="${STAGE}/per_sample/${id}/${id}_stats.txt"
        require_nonempty_file "${flag}"
        require_nonempty_file "${stats}"

        total="$(awk '/ in total / {print $1; exit}' "${flag}")"
        mapped_line="$(grep -m 1 ' mapped (' "${flag}")"
        proper_line="$(grep -m 1 ' properly paired (' "${flag}")"
        mapped="${mapped_line%% *}"
        proper="${proper_line%% *}"
        mapped_pct="$(sed -n 's/.*(\([0-9.][0-9.]*\)%.*/\1/p' <<< "${mapped_line}")"
        proper_pct="$(sed -n 's/.*(\([0-9.][0-9.]*\)%.*/\1/p' <<< "${proper_line}")"
        duplicates="$(awk '/ duplicates$/ {print $1; exit}' "${flag}")"
        duplicate_pct="$(awk -v d="${duplicates:-0}" -v t="${total:-0}" 'BEGIN {if (t>0) printf "%.4f", 100*d/t; else print ""}')"
        bases_cigar="$(awk -F '\t' '$1=="SN" && $2 ~ /^bases mapped \(cigar\):/ {print $3; exit}' "${stats}")"

        printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
            "${id}" "${total}" "${mapped}" "${mapped_pct}" \
            "${proper}" "${proper_pct}" "${duplicates}" "${duplicate_pct}" "${bases_cigar}" \
            >> "${summary_tmp}"
    done < "${MANIFEST}"

    mv "${summary_tmp}" "${STAGE}/summary_metrics.csv"
    require_nonempty_file "${STAGE}/summary_metrics.csv"

    mv "${STAGE}" "${ROOT}/05_quality"
fi

# ==============================================================================
# 11_featureCounts — paired-end gene counts with Rsubread 2.16.1
# This is the supplied draft endpoint. Sample columns are explicitly renamed
# from the manifest IDs, replacing the draft's fragile path-prefix/suffix logic.
# ==============================================================================
if begin_stage "11_featureCounts"; then
    echo "Creating and running 11_featureCounts"
    STAGE="$(mktemp -d "${ROOT}/11_featureCounts.tmp.XXXXXX")"

    ensure_sif \
        "rsubread_2.16.1.sif" \
        "https://depot.galaxyproject.org/singularity/bioconductor-rsubread%3A2.16.1--r43ha9d7317_0" \
        "docker://quay.io/biocontainers/bioconductor-rsubread:2.16.1--r43ha9d7317_0"

    bam_manifest_tmp="$(mktemp "${STAGE}/.bam_manifest.tsv.tmp.XXXXXX")"
    printf 'ID\tBAM\n' > "${bam_manifest_tmp}"

    while IFS=, read -r id r1 r2; do
        id="${id//$'\r'/}"
        [[ "${id}" == "ID" || -z "${id}" ]] && continue
        bam="${ROOT}/05_MarkDuplicates/bam/${id}.mkdup.sorted.bam"
        require_nonempty_file "${bam}"
        printf '%s\t%s\n' "${id}" "${bam}" >> "${bam_manifest_tmp}"
    done < "${MANIFEST}"

    mv "${bam_manifest_tmp}" "${STAGE}/bam_manifest.tsv"

    cat > "${STAGE}/featureCounts.R" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
    stop("Usage: featureCounts.R <annotation.gtf> <bam_manifest.tsv> <cnt.csv> <assignment_summary.csv>")
}

annotation <- args[[1L]]
bam_manifest <- args[[2L]]
counts_out <- args[[3L]]
summary_out <- args[[4L]]

suppressPackageStartupMessages(library(Rsubread))

if (as.character(packageVersion("Rsubread")) != "2.16.1") {
    stop("Rsubread version mismatch: expected 2.16.1, found ", packageVersion("Rsubread"))
}

samples <- read.delim(
    bam_manifest,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (nrow(samples) == 0L || anyDuplicated(samples$ID) || any(!file.exists(samples$BAM))) {
    stop("Invalid BAM manifest")
}

fc <- featureCounts(
    files = samples$BAM,
    annot.ext = annotation,
    isGTFAnnotationFile = TRUE,
    isPairedEnd = TRUE,
    nthreads = 4
)

colnames(fc$counts) <- samples$ID

count_table <- data.frame(
    fc$annotation[, c("GeneID", "Length"), drop = FALSE],
    fc$counts,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

write.table(
    count_table,
    file = counts_out,
    quote = FALSE,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE
)

assignment <- as.data.frame(fc$stat, check.names = FALSE, stringsAsFactors = FALSE)
colnames(assignment)[-1L] <- samples$ID
write.table(
    assignment,
    file = summary_out,
    quote = FALSE,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE
)

writeLines(
    c(
        paste0("Rsubread=", packageVersion("Rsubread")),
        capture.output(sessionInfo())
    ),
    con = file.path(dirname(counts_out), "sessionInfo.txt")
)
RSCRIPT

    submission="$(
        ROOT="${ROOT}" STAGE="${STAGE}" ANNOTATION="${ANNOTATION}" RUNTIME="${CONTAINER_RUNTIME}" RSUBREAD_SIF="${RSUBREAD_SIF}" \
        bsub \
            -J "featureCounts" \
            -n 4 \
            -R "select[mem>=16000] rusage[mem=16000] span[hosts=1]" \
            -M 16000 \
            -o "${STAGE}/featureCounts.%J.out" \
            -e "${STAGE}/featureCounts.%J.err" \
            bash -lc '
set -euo pipefail

counts_tmp="$(mktemp "${STAGE}/.cnt.csv.tmp.XXXXXX")"
summary_tmp="$(mktemp "${STAGE}/.featureCounts_assignment_summary.csv.tmp.XXXXXX")"

"${RUNTIME}" exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${RSUBREAD_SIF}" \
    Rscript \
    "${STAGE}/featureCounts.R" \
    "${ANNOTATION}" \
    "${STAGE}/bam_manifest.tsv" \
    "${counts_tmp}" \
    "${summary_tmp}"

[[ -s "${counts_tmp}" && -s "${summary_tmp}" ]]

expected_cols="$(( $(awk '"'"'END {print NR-1}'"'"' "${STAGE}/bam_manifest.tsv") + 2 ))"
observed_cols="$(awk -F',' '"'"'NR==1 {print NF; exit}'"'"' "${counts_tmp}")"
[[ "${observed_cols}" -eq "${expected_cols}" ]]

mv "${counts_tmp}" "${STAGE}/cnt.csv"
mv "${summary_tmp}" "${STAGE}/featureCounts_assignment_summary.csv"
'
    )"

    echo "${submission}"
    jid="$(extract_job_id "${submission}")"
    wait_for_job "${jid}"

    require_nonempty_file "${STAGE}/cnt.csv"
    require_nonempty_file "${STAGE}/featureCounts_assignment_summary.csv"

    mv "${STAGE}" "${ROOT}/11_featureCounts"
fi

echo "Transcriptomics pipeline completed successfully."
