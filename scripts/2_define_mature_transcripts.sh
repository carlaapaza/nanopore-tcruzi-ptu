#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Script: 2_define_mature_transcripts.sh
# Purpose: Identify SL-insertion (5') and polyA-addition (3') sites from
#          ONT direct RNA-seq (stranded) in Trypanosomatids (e.g., T. cruzi).
#
# Repo layout this script expects (run from repo ROOT):
#   data/genome/yc6.fna
#   data/reads/Tcruzi_ystrain.fastq.gz
#   data/annotation/ptu.gff           (optional; set PTU_GFF below)
#   results/02_mature_transcripts/    (outputs created here)
#
# Pipeline
#   1) Detect reads containing the 39-nt Spliced Leader (SL) at 5' (anchored).
#   2) Map SL-trimmed reads to genome (minimap2) → sorted/indexed BAM.
#   3) Call strand-aware SL insertion sites (BED) + counts (TSV).
#   4) Detect reads with a 3' polyA tail (anchored), trim, map to genome.
#   5) Call strand-aware polyA addition sites (BED) + counts (TSV).
#   6) (Optional) Check SL orientation vs PTUs (if PTU_GFF is provided).
#
#
# Requirements:
#   - cutadapt, minimap2, samtools, bedtools, (seqkit optional)
#
# Author: Carla Apaza
# ============================================================

# ------------------------ CONFIG (edit paths or export as env) ----------------------------
# Inputs (relative to repo root)
READS="${READS:-data/reads/Tcruzi_ystrain.fastq.gz}"     # ONT direct RNA FASTQ(.gz)
GENOME="${GENOME:-data/genome/yc6.fna}"                  # Genome FASTA
PTU_GFF="${PTU_GFF:-}"                                   # Optional: data/annotation/ptu.gff

# Output folder
OUTDIR="${OUTDIR:-results/02_mature_transcripts}"
TMPDIR="${OUTDIR}/tmp"
mkdir -p "${OUTDIR}" "${TMPDIR}" "${OUTDIR}/logs"

# SL sequence (39 nt) — your T. cruzi SL
SL_SEQ="${SL_SEQ:-AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG}"

# Detection thresholds
SL_MIN_OVERLAP="${SL_MIN_OVERLAP:-20}"     # min overlap with SL motif
SL_MAX_ERROR="${SL_MAX_ERROR:-0.22}"       # max error rate (ONT-tolerant)
POLYA_MIN_RUN="${POLYA_MIN_RUN:-15}"       # min A's to call a 3' polyA tail

# Mapping preset / threads
MAP_PRESET="${MAP_PRESET:-map-ont}"
THREADS="${THREADS:-12}"

# Conda environments (edit to your env names)
CUTADAPT_ENV="${CUTADAPT_ENV:-cutadapt}"
MINIMAP2_ENV="${MINIMAP2_ENV:-minimap2}"
SAMTOOLS_ENV="${SAMTOOLS_ENV:-samtools}"
BEDTOOLS_ENV="${BEDTOOLS_ENV:-bedtools_env}"
SEQKIT_ENV="${SEQKIT_ENV:-seqkit}"

# ------------------------ Logging helpers -------------------
msg()  { echo -e "[$(date +'%F %T')] $*"; }
fail() { echo -e "[$(date +'%F %T')][FATAL] $*" >&2; exit 1; }

# ------------------------ Conda init ------------------------
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh" || fail "conda not found in PATH"

activate() {
  local env_name="$1"
  msg "Activating conda env: ${env_name}"
  conda activate "${env_name}" || fail "conda activate ${env_name} failed"
}

# ------------------------ Guards ----------------------------
[[ -s "${GENOME}" ]] || fail "Genome not found: ${GENOME}"
[[ -s "${READS}"  ]] || fail "Reads not found: ${READS}"

# ------------------------ Idempotent helpers ----------------
run_if_missing() {
  local outfile="$1"; shift
  local env="$1"; shift
  if [[ -s "${outfile}" ]]; then
    msg "[SKIP] ${outfile}"
  else
    msg "[RUN ] Producing ${outfile}"
    activate "${env}"
    # shellcheck disable=SC2068
    "$@" || fail "Command failed for ${outfile}"
    [[ -s "${outfile}" ]] || fail "${outfile} not created or empty."
  fi
}

run_block() {
  local marker="$1"; shift
  local env="$1"; shift
  local desc="$1"; shift
  if [[ -s "${marker}" ]]; then
    msg "[SKIP] ${desc} (marker: $(basename "${marker}"))"
  else
    msg "[RUN ] ${desc}"
    activate "${env}"
    local __code; __code="$(cat)"
    eval "${__code}" || fail "Block failed: ${desc}"
    : > "${marker}"
  fi
}

# ------------------------ PREP: minimap2 index --------------
GENOME_MMI="${GENOME}.mmi"
if [[ ! -s "${GENOME_MMI}" ]]; then
  msg "Building minimap2 index for genome..."
  activate "${MINIMAP2_ENV}"
  minimap2 -d "${GENOME_MMI}" "${GENOME}"
else
  msg "[SKIP] Genome index exists: ${GENOME_MMI}"
fi

# ------------------------ SL fasta (forward-only; stranded) -
SL_FA="${OUTDIR}/sl_5prime.fa"
if [[ ! -s "${SL_FA}" ]]; then
  msg "Creating SL fasta (forward-only): ${SL_FA}"
  { echo ">SL"; echo "${SL_SEQ}"; } > "${SL_FA}"
else
  msg "[SKIP] SL fasta exists: ${SL_FA}"
fi

# ------------------------ STEP 1: 5'-anchored SL (cutadapt) -
SL_TRIM_FQ="${OUTDIR}/reads_with_SL.trimmed.fastq.gz"
SL_IDS="${OUTDIR}/reads_with_SL.ids"
run_block "${OUTDIR}/.step1.done" "${CUTADAPT_ENV}" "Detecting 5'-anchored SL with cutadapt (stranded)" <<EOF
cutadapt \
  -g "^file:${SL_FA}" \
  --overlap ${SL_MIN_OVERLAP} \
  -e ${SL_MAX_ERROR} \
  --discard-untrimmed \
  -o "${SL_TRIM_FQ}" \
  "${READS}" \
  > "${OUTDIR}/logs/cutadapt_SL.log" 2>&1

# Extract read IDs (those that had SL trimmed)
if [[ "${SL_TRIM_FQ}" == *.gz ]]; then
  zcat "${SL_TRIM_FQ}" | awk 'NR%4==1 {sub(/^@/,""); print \$1}' > "${SL_IDS}"
else
  awk 'NR%4==1 {sub(/^@/,""); print \$1}' "${SL_TRIM_FQ}" > "${SL_IDS}"
fi
EOF

if [[ ! -s "${SL_TRIM_FQ}" ]]; then
  msg "[WARN] No SL-containing reads found. Check SL sequence/parameters."
fi

# ------------------------ STEP 2: map SL-trimmed → BAM ------
SL_SAM="${OUTDIR}/sl.trimmed.toGenome.sam"
SL_BAM="${OUTDIR}/sl.trimmed.toGenome.bam"
SL_SORT_BAM="${OUTDIR}/sl.trimmed.toGenome.sorted.bam"
SL_SORT_BAI="${SL_SORT_BAM}.bai"

# 2a) minimap2 → SAM
if [[ -s "${SL_TRIM_FQ}" ]]; then
  run_if_missing "${SL_SAM}" "${MINIMAP2_ENV}" \
    minimap2 -t "${THREADS}" -ax "${MAP_PRESET}" "${GENOME_MMI}" "${SL_TRIM_FQ}" -o "${SL_SAM}"
else
  msg "[WARN] ${SL_TRIM_FQ} empty; skipping mapping."
  : > "${SL_SAM}"
fi

# 2b) SAM → BAM
if [[ -s "${SL_SAM}" ]]; then
  run_if_missing "${SL_BAM}" "${SAMTOOLS_ENV}" \
    samtools view -@ "${THREADS}" -b -F 4 -o "${SL_BAM}" "${SL_SAM}"
fi

# 2c) sort + index
if [[ -s "${SL_BAM}" ]]; then
  run_if_missing "${SL_SORT_BAM}" "${SAMTOOLS_ENV}" \
    samtools sort -@ "${THREADS}" -o "${SL_SORT_BAM}" "${SL_BAM}"
  if [[ ! -s "${SL_SORT_BAI}" ]]; then
    msg "[RUN ] Indexing SL BAM ..."
    activate "${SAMTOOLS_ENV}"
    samtools index "${SL_SORT_BAM}"
  else
    msg "[SKIP] ${SL_SORT_BAI} exists."
  fi
fi

# ------------------------ STEP 3: SL sites -------------------
SL_BED_RAW="${OUTDIR}/sl.alignments.bed"
SL_SITES_BED="${OUTDIR}/SL_insertion_sites.bed"
SL_SITES_TSV="${OUTDIR}/SL_insertion_sites.counts.tsv"

run_if_missing "${SL_BED_RAW}" "${BEDTOOLS_ENV}" \
  bedtools bamtobed -i "${SL_SORT_BAM}" > "${SL_BED_RAW}"

run_block "${OUTDIR}/.step3.done" "${BEDTOOLS_ENV}" "Computing SL insertion sites (5' ends)" <<'EOF'
# 1-bp site at the 5' end of each alignment (strand-aware)
awk 'BEGIN{OFS="\t"}{
  chrom=$1; start=$2; end=$3; strand=$6;
  if (strand=="+") { s=start; e=start+1 }
  else if (strand=="-") { s=end-1; e=end }
  else { next }
  print chrom, s, e, ".", 0, strand
}' "${SL_BED_RAW}" > "${SL_SITES_BED}"

# Collapse to unique sites and count
awk 'BEGIN{OFS="\t"}{key=$1 FS $2 FS $3 FS $6; c[key]++}END{
  print "chrom","start","end","count","strand";
  for(k in c){split(k,a,FS); print a[1],a[2],a[3],c[k],a[4]}
}' "${SL_SITES_BED}" \
| sort -k1,1 -k2,2n > "${SL_SITES_TSV}"
EOF

# ------------------------ STEP 4: 3'-anchored polyA ----------
POLYA_TRIM_FQ="${OUTDIR}/reads_with_polyA.trimmed.fastq.gz"
POLYA_SAM="${OUTDIR}/polyA.trimmed.toGenome.sam"
POLYA_BAM="${OUTDIR}/polyA.trimmed.toGenome.bam"
POLYA_SORT_BAM="${OUTDIR}/polyA.trimmed.toGenome.sorted.bam"
POLYA_SORT_BAI="${POLYA_SORT_BAM}.bai"

# 4a) cutadapt (anchor tail)
run_if_missing "${POLYA_TRIM_FQ}" "${CUTADAPT_ENV}" \
  bash -lc "cutadapt -a \"A{${POLYA_MIN_RUN}}$\" --discard-untrimmed -o \"${POLYA_TRIM_FQ}\" \"${READS}\" > \"${OUTDIR}/logs/cutadapt_polyA.log\" 2>&1"

# 4b) minimap2 → SAM
run_if_missing "${POLYA_SAM}" "${MINIMAP2_ENV}" \
  minimap2 -t "${THREADS}" -ax "${MAP_PRESET}" "${GENOME_MMI}" "${POLYA_TRIM_FQ}" -o "${POLYA_SAM}"

# 4c) SAM → BAM
run_if_missing "${POLYA_BAM}" "${SAMTOOLS_ENV}" \
  samtools view -@ "${THREADS}" -b -F 4 -o "${POLYA_BAM}" "${POLYA_SAM}"

# 4d) sort + index
run_if_missing "${POLYA_SORT_BAM}" "${SAMTOOLS_ENV}" \
  samtools sort -@ "${THREADS}" -o "${POLYA_SORT_BAM}" "${POLYA_BAM}"

if [[ -s "${POLYA_SORT_BAM}" && ! -s "${POLYA_SORT_BAI}" ]]; then
  msg "[RUN ] Indexing polyA BAM ..."
  activate "${SAMTOOLS_ENV}"
  samtools index "${POLYA_SORT_BAM}"
fi

# ------------------------ STEP 5: polyA sites ----------------
POLYA_BED_RAW="${OUTDIR}/polyA.alignments.bed"
POLYA_SITES_BED="${OUTDIR}/polyA_addition_sites.bed"
POLYA_SITES_TSV="${OUTDIR}/polyA_addition_sites.counts.tsv"

run_if_missing "${POLYA_BED_RAW}" "${BEDTOOLS_ENV}" \
  bedtools bamtobed -i "${POLYA_SORT_BAM}" > "${POLYA_BED_RAW}"

run_block "${OUTDIR}/.step5.done" "${BEDTOOLS_ENV}" "Computing polyA addition sites (3' ends)" <<'EOF'
# 1-bp site at the 3' end (strand-aware)
awk 'BEGIN{OFS="\t"}{
  chrom=$1; start=$2; end=$3; strand=$6;
  if (strand=="+") { s=end-1; e=end }
  else if (strand=="-") { s=start; e=start+1 }
  else { next }
  print chrom, s, e, ".", 0, strand
}' "${POLYA_BED_RAW}" > "${POLYA_SITES_BED}"

# Collapse to unique sites and count
awk 'BEGIN{OFS="\t"}{key=$1 FS $2 FS $3 FS $6; c[key]++}END{
  print "chrom","start","end","count","strand";
  for(k in c){split(k,a,FS); print a[1],a[2],a[3],c[k],a[4]}
}' "${POLYA_SITES_BED}" \
| sort -k1,1 -k2,2n > "${POLYA_SITES_TSV}"
EOF

# ------------------------ STEP 6: Summaries ------------------
msg "Summaries"
activate "${SAMTOOLS_ENV}"
TOTAL_READS=$( ( [[ "${READS}" == *.gz ]] && zcat "${READS}" || cat "${READS}" ) | awk 'END{print NR/4}' )
SL_READS=$( ( [[ "${SL_TRIM_FQ}" == *.gz ]] && zcat "${SL_TRIM_FQ}" || cat "${SL_TRIM_FQ}" ) | awk 'END{print NR/4}' )
POLYA_READS=$( ( [[ "${POLYA_TRIM_FQ}" == *.gz ]] && zcat "${POLYA_TRIM_FQ}" || cat "${POLYA_TRIM_FQ}" ) | awk 'END{print NR/4}' )
SL_SITES=$(awk 'NR>1' "${SL_SITES_TSV}" | wc -l | awk '{print $1}')
POLYA_SITES=$(awk 'NR>1' "${POLYA_SITES_TSV}" | wc -l | awk '{print $1}')

cat <<SUM
=== Summary ===
  Total reads:                      ${TOTAL_READS}
  Reads with SL motif (5'-anchored): ${SL_READS}
  Unique SL insertion sites:        ${SL_SITES}
  Reads with polyA tail (3'-anch.): ${POLYA_READS}
  Unique polyA addition sites:      ${POLYA_SITES}
Logs:
  - ${OUTDIR}/logs/cutadapt_SL.log
  - ${OUTDIR}/logs/cutadapt_polyA.log
SUM

# ------------------------ STEP 7: SL vs PTUs (optional) -----
if [[ -n "${PTU_GFF}" && -s "${PTU_GFF}" ]]; then
  msg "Checking SL orientation vs PTUs ..."
  PTU_BED="${OUTDIR}/ptu.bed"
  activate "${BEDTOOLS_ENV}"
  awk 'BEGIN{OFS="\t"} $0 !~ /^#/ {
      chrom=$1; feat=$3; start=$4-1; end=$5; strand=$7; name=(NF>=9?$9:"PTU");
      if (feat=="PTU" || feat=="transcription_unit" || feat=="gene_cluster") print chrom,start,end,name,0,strand
  }' "${PTU_GFF}" | sort -k1,1 -k2,2n > "${PTU_BED}"

  SL_IN_PTU="${OUTDIR}/SL_sites.inPTU.tsv"
  bedtools intersect -s -wa -wb -a "${SL_SITES_BED}" -b "${PTU_BED}" \
    | awk 'BEGIN{OFS="\t"; print "chrom","start","end","SL_strand","PTU_name","PTU_strand"}{print $1,$2,$3,$6,$10,$12}' > "${SL_IN_PTU}"

  CONTRARY="${OUTDIR}/SL_sites.contrary_orientation.tsv"
  bedtools intersect -wa -wb -a "${SL_SITES_BED}" -b "${PTU_BED}" \
    | awk 'BEGIN{OFS="\t"; print "chrom","start","end","SL_strand","PTU_name","PTU_strand"}{ if ($6!=$12) print $1,$2,$3,$6,$10,$12 }' > "${CONTRARY}"

  msg "SL sites within PTUs (same strand):  $(($(wc -l < "${SL_IN_PTU}")-1))"
  msg "SL sites contrary within PTUs:       $(($(wc -l < "${CONTRARY}")-1))  (should be ~0)"
else
  msg "PTU_GFF not provided — skipping orientation check."
fi

msg "[DONE] Outputs in: ${OUTDIR}"
msg "  - ${SL_SITES_TSV}  (SL insertion sites with counts)"
msg "  - ${POLYA_SITES_TSV} (polyA addition sites with counts)"
msg "  - ${SL_SITES_BED} , ${POLYA_SITES_BED} (1-bp BED, strand-aware)"
