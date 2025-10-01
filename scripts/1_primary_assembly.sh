#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Script: 1_primary_assembly.sh
# Purpose: Genome-guided PRIMARY transcriptome assembly
#          for ONT direct RNA-seq (T. cruzi Y clone 6).
#
# Produces a "primary assembly" of polycistronic units (PTUs)
# to be later refined with SL/polyA site evidence.
#
# Pipeline:
#   minimap2 (spliced) → sorted BAM → StringTie2 (conservative)
#   → gffread FASTAs
#
# Author: Carla Apaza (T. cruzi project)
# ============================================================

# --------------------------- CONFIG ---------------------------
REF_GENOME="data/genome/yc6.fna"
READS_FASTQ="data/reads/Tcruzi_ystrain.fastq.gz"
OUTDIR="results/1_primary_assembly"
THREADS=8

REF_ANNOT_GTF=""

# Conda envs
ENV_MINIMAP2="minimap2"
ENV_SAMTOOLS="samtools"
ENV_STRINGTIE="stringtie"
ENV_GFFREAD="gffread"
ENV_SEQKIT="seqkit"

# StringTie conservative params (PTU-like)
STRINGTIE_MINLEN=200
STRINGTIE_F=0.15
STRINGTIE_J=2
STRINGTIE_C=2
STRINGTIE_EXTRA="--conservative"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# --------------------------- Logging ---------------------------
msg() { echo -e "[$(date +'%F %T')] $*"; }
fail() { echo -e "[$(date +'%F %T')][FATAL] $*" >&2; exit 1; }

# --------------------------- Conda init ---------------------------
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

# --------------------------- Checks ---------------------------
[[ -s "${REF_GENOME}" ]] || fail "Genome not found: ${REF_GENOME}"
[[ -s "${READS_FASTQ}" ]] || fail "Reads not found: ${READS_FASTQ}"

# --------------------------- Paths ---------------------------
REF_IDX="${OUTDIR}/genome_ref.mmi"
ALN_SAM="${OUTDIR}/aln.genome.spliced.sam"
ALN_BAM="${OUTDIR}/aln.genome.spliced.bam"
ALN_BAM_SORTED="${OUTDIR}/aln.genome.spliced.sorted.bam"
ALN_BAI="${ALN_BAM_SORTED}.bai"
FLAGSTAT="${OUTDIR}/aln.flagstat.txt"

STRG_GTF="${OUTDIR}/stringtie.primary_ptu.gtf"
STRG_ABUND_TSV="${OUTDIR}/stringtie.primary_ptu.abund.tsv"
TX_FASTA="${OUTDIR}/transcripts.primary_ptu.fa"
CDS_FASTA="${OUTDIR}/transcripts.primary_ptu.cds.fa"

msg "Genome : ${REF_GENOME}"
msg "Reads  : ${READS_FASTQ}"
msg "Outdir : ${OUTDIR}"
msg "StringTie (PTU-like): -L ${STRINGTIE_EXTRA} -m ${STRINGTIE_MINLEN} -f ${STRINGTIE_F} -j ${STRINGTIE_J} -c ${STRINGTIE_C}"

# --------------------------- 1) minimap2 index ---------------------------
if [[ -s "${REF_IDX}" ]]; then
  msg "[SKIP] Index exists: ${REF_IDX}"
else
  msg "[RUN] Building minimap2 index..."
  conda activate "${ENV_MINIMAP2}"
  minimap2 -d "${REF_IDX}" "${REF_GENOME}"
  conda deactivate || true
fi

# --------------------------- 2) spliced alignment ---------------------------
if [[ -s "${ALN_SAM}" ]]; then
  msg "[SKIP] SAM exists: ${ALN_SAM}"
else
  msg "[RUN] Spliced alignment (minimap2 -> SAM)..."
  conda activate "${ENV_MINIMAP2}"
  minimap2 -t "${THREADS}" -ax splice -uf -k14 "${REF_IDX}" "${READS_FASTQ}" > "${ALN_SAM}"
  conda deactivate || true
fi

# --------------------------- 3) SAM/BAM conversion ---------------------------
if [[ -s "${ALN_BAM_SORTED}" && -s "${ALN_BAI}" ]]; then
  msg "[SKIP] Sorted BAM and index exist: ${ALN_BAM_SORTED}"
else
  msg "[RUN] Converting, sorting and indexing (samtools)..."
  conda activate "${ENV_SAMTOOLS}"

  if [[ -s "${ALN_BAM}" ]]; then
    msg "[SKIP] Unsorted BAM exists: ${ALN_BAM}"
  else
    samtools view -@ "${THREADS}" -b "${ALN_SAM}" -o "${ALN_BAM}"
  fi

  samtools sort  -@ "${THREADS}" -o "${ALN_BAM_SORTED}" "${ALN_BAM}"
  samtools index "${ALN_BAM_SORTED}"
  conda deactivate || true
fi

# QC metrics
if [[ -s "${FLAGSTAT}" ]]; then
  msg "[SKIP] Flagstat exists: ${FLAGSTAT}"
else
  msg "[RUN] Computing flagstat..."
  conda activate "${ENV_SAMTOOLS}"
  samtools flagstat "${ALN_BAM_SORTED}" > "${FLAGSTAT}"
  conda deactivate || true
fi

# --------------------------- 4) StringTie primary assembly ---------------------------
if [[ -s "${STRG_GTF}" ]]; then
  msg "[SKIP] StringTie GTF exists: ${STRG_GTF}"
else
  msg "[RUN] StringTie2 primary assembly..."
  conda activate "${ENV_STRINGTIE}"
  stringtie -L -p "${THREADS}" ${STRINGTIE_EXTRA} \
    -m "${STRINGTIE_MINLEN}" -f "${STRINGTIE_F}" -j "${STRINGTIE_J}" -c "${STRINGTIE_C}" \
    -A "${STRG_ABUND_TSV}" \
    -o "${STRG_GTF}" \
    "${ALN_BAM_SORTED}"
  conda deactivate || true
fi

# --------------------------- 5) export FASTAs ---------------------------
if [[ -s "${TX_FASTA}" ]]; then
  msg "[SKIP] Transcript FASTA exists: ${TX_FASTA}"
else
  msg "[RUN] Exporting transcript FASTA..."
  conda activate "${ENV_GFFREAD}"
  gffread -w "${TX_FASTA}" -g "${REF_GENOME}" "${STRG_GTF}"
  conda deactivate || true
fi

if [[ -s "${CDS_FASTA}" ]]; then
  msg "[SKIP] CDS FASTA exists: ${CDS_FASTA}"
else
  msg "[RUN] Exporting CDS FASTA..."
  conda activate "${ENV_GFFREAD}"
  gffread -x "${CDS_FASTA}" -g "${REF_GENOME}" "${STRG_GTF}" || true
  conda deactivate || true
fi

msg "[DONE] PRIMARY PTU-like assembly complete."
msg " - Sorted BAM/index : ${ALN_BAM_SORTED} , ${ALN_BAI}"
msg " - Mapping metrics  : ${FLAGSTAT}"
msg " - Primary GTF      : ${STRG_GTF}"
msg " - Abundance (TSV)  : ${STRG_ABUND_TSV}"
msg " - Transcripts FASTA: ${TX_FASTA}"
msg " - CDS FASTA (opt)  : ${CDS_FASTA}"
