#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Script: genome_guided_transcriptome.sh
# Purpose: Genome-guided transcriptome assembly for ONT direct RNA
# Pipeline:
#   minimap2 (spliced) -> BAM -> StringTie2 (-L, m,f,j,c) -> gffread -> transcripts.fa (+CDS)
# Notes:
#   - Idempotent: steps are skipped if outputs already exist.
#   - Uses conda env switching with direct `conda activate` / `conda deactivate`.
# Author: Carla Apaza (T. cruzi project)
# ============================================================

# --------------------------- Logging ---------------------------
msg() { echo -e "[$(date +'%F %T')] $*"; }

# ------------------------ Conda init --------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
export QT_XCB_GL_INTEGRATION="${QT_XCB_GL_INTEGRATION:-}"

# -------------------------- Defaults --------------------------
REF_GENOME="/home/block/Desktop/02-PTUTcruzi/final/ref_assembly_genome/yc6/yc6.fna"
READS_FASTQ="/home/block/Desktop/02-PTUTcruzi/data/Tcruzi_ystrain.fastq.gz"
OUTDIR="/home/block/Desktop/02-PTUTcruzi/final/ref_assembly_genome/assembly_40"
THREADS=8
REF_ANNOT_GTF="/home/block/Desktop/02-PTUTcruzi/final/ref_assembly_genome/yc6/genomic.gff"

# Conda environments
ENV_MINIMAP2="minimap2"
ENV_SAMTOOLS="samtools"
ENV_STRINGTIE="stringtie"
ENV_GFFREAD="gffread"
ENV_SEQKIT="seqkit"

# StringTie params
STRINGTIE_MINLEN=40
STRINGTIE_F=0.01
STRINGTIE_J=1
STRINGTIE_C=1

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# -------------------------- Paths -----------------------------
REF_IDX="${OUTDIR}/genome_ref.mmi"
ALN_BAM_SORTED="${OUTDIR}/aln.genome.spliced.sorted.bam"
ALN_BAI="${ALN_BAM_SORTED}.bai"
STRG_GTF="${OUTDIR}/stringtie.transcripts.gtf"
STRG_ABUND_TSV="${OUTDIR}/stringtie.gene_abund.tsv"
TX_FASTA="${OUTDIR}/transcripts.genome_guided.fa"
CDS_FASTA="${OUTDIR}/transcripts.genome_guided.cds.fa"

# =================== 1) minimap2 index =======================
if [[ -s "${REF_IDX}" ]]; then
  msg "[SKIP] Index exists: ${REF_IDX}"
else
  msg "[RUN ] Building minimap2 index..."
  conda activate "${ENV_MINIMAP2}"
  minimap2 -d "${REF_IDX}" "${REF_GENOME}"
  conda deactivate
fi

# ===== 2) Spliced alignment (ONT) -> sorted BAM ==============
if [[ -s "${ALN_BAM_SORTED}" && -s "${ALN_BAI}" ]]; then
  msg "[SKIP] Sorted BAM exists: ${ALN_BAM_SORTED}"
else
  msg "[RUN ] Spliced alignment with minimap2..."
  conda activate "${ENV_MINIMAP2}"
  minimap2 -t "${THREADS}" -ax splice -uf -k14 "${REF_IDX}" "${READS_FASTQ}" \
    | samtools sort -@ "${THREADS}" -o "${ALN_BAM_SORTED}"
  conda deactivate

  conda activate "${ENV_SAMTOOLS}"
  samtools index "${ALN_BAM_SORTED}"
  conda deactivate
fi

# === 3) Reference-guided assembly (StringTie2) ===============
if [[ -s "${STRG_GTF}" ]]; then
  msg "[SKIP] StringTie GTF exists: ${STRG_GTF}"
else
  msg "[RUN ] StringTie2 assembly..."
  conda activate "${ENV_STRINGTIE}"
  if [[ -n "${REF_ANNOT_GTF}" ]]; then
    stringtie -L -p "${THREADS}" \
      -m "${STRINGTIE_MINLEN}" -f "${STRINGTIE_F}" -j "${STRINGTIE_J}" -c "${STRINGTIE_C}" \
      -G "${REF_ANNOT_GTF}" \
      -A "${STRG_ABUND_TSV}" \
      -o "${STRG_GTF}" \
      "${ALN_BAM_SORTED}"
  else
    stringtie -L -p "${THREADS}" \
      -m "${STRINGTIE_MINLEN}" -f "${STRINGTIE_F}" -j "${STRINGTIE_J}" -c "${STRINGTIE_C}" \
      -A "${STRG_ABUND_TSV}" \
      -o "${STRG_GTF}" \
      "${ALN_BAM_SORTED}"
  fi
  conda deactivate
fi

# =========== 4) Export transcript FASTA and CDS ==============
if [[ -s "${TX_FASTA}" ]]; then
  msg "[SKIP] Transcript FASTA exists: ${TX_FASTA}"
else
  msg "[RUN ] Exporting transcript FASTA with gffread..."
  conda activate "${ENV_GFFREAD}"
  gffread -w "${TX_FASTA}" -g "${REF_GENOME}" "${STRG_GTF}"
  conda deactivate
fi

if [[ -s "${CDS_FASTA}" ]]; then
  msg "[SKIP] CDS FASTA exists: ${CDS_FASTA}"
else
  msg "[RUN ] Exporting CDS FASTA with gffread..."
  conda activate "${ENV_GFFREAD}"
  gffread -x "${CDS_FASTA}" -g "${REF_GENOME}" "${STRG_GTF}" || true
  conda deactivate
fi

# =========== 5) FASTA stats (optional) =======================
if [[ -s "${OUTDIR}/transcripts.stats.txt" ]]; then
  msg "[SKIP] FASTA stats exist."
else
  if command -v seqkit >/dev/null 2>&1; then
    msg "[RUN ] Computing FASTA stats with seqkit..."
    conda activate "${ENV_SEQKIT}"
    seqkit stats -a "${TX_FASTA}" > "${OUTDIR}/transcripts.stats.txt" || true
    conda deactivate
  else
    msg "[INFO] Skipping seqkit stats (seqkit not available)."
  fi
fi

# -------------------------- Summary --------------------------
msg "[DONE] Genome-guided transcriptome assembly completed."
msg " - Sorted BAM/index : ${ALN_BAM_SORTED} , ${ALN_BAI}"
msg " - Transcripts (GTF): ${STRG_GTF}"
msg " - Abundance (TSV)  : ${STRG_ABUND_TSV}"
msg " - Transcripts FASTA: ${TX_FASTA}"
msg " - CDS FASTA (opt)  : ${CDS_FASTA}"
