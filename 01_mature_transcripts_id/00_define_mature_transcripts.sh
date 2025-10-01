#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Script: define_mature_transcripts.sh
# Purpose: Identify SL-insertion (5') and polyA-addition (3') sites from
#          ONT direct RNA-seq in Trypanosomatids (e.g., T. cruzi).
#
# Pipeline (high-level):
#   1) Detect reads containing the 39-nt Spliced Leader (SL) motif anywhere.
#      - Trim the SL motif (ONT-error tolerant).
#   2) Map SL-trimmed reads to the genome (minimap2) -> sorted/indexed BAM.
#   3) Call strand-aware SL insertion sites (1-bp BED) + counts (TSV).
#   4) Detect reads with 3' polyA tails (>=N As), trim, map to genome.
#   5) Call strand-aware polyA addition sites (1-bp BED) + counts (TSV).
#   6) (Optional) Check SL orientation vs PTU annotation (GFF/GTF).
#
# Idempotency:
#   - Each step skips if the expected output already exists and is non-empty.
#     Delete an output (or the step marker) to force re-run.
#
# Conda:
#   - Before each tool, activates its conda env (configurable below).
#
# Requirements (in their respective conda envs):
#   - cutadapt >=4, minimap2, samtools, bedtools
#   - seqkit (optional; not required)
#
# Author: Carla Apaza (T. cruzi project)
# ============================================================

# --------------------------- Logging ---------------------------
msg() { echo -e "[$(date +'%F %T')] $*"; }
fail() { echo -e "[FATAL] $*" >&2; exit 1; }
trap 'echo "[FATAL] Line $LINENO (exit=$?)" >&2' ERR

# --------------------------- Usage ----------------------------
usage() {
  cat <<EOF
Usage:
  $(basename "$0") [--reads FASTQ(.gz)] [--genome FASTA] [--outdir DIR]
                   [--ptu-gff FILE] [--threads N]
                   [--sl-seq SEQ39] [--sl-min-overlap N] [--sl-max-error F]
                   [--polya-min-run N]
                   [--cutadapt-env NAME] [--minimap2-env NAME]
                   [--samtools-env NAME] [--bedtools-env NAME]

Example:
  $(basename "$0") \\
    --reads /path/Tcruzi_ystrain.fastq.gz \\
    --genome /path/yc6.fna \\
    --outdir /path/02_mature_transcripts \\
    --ptu-gff /path/ptus.gff3

Notes:
  - Outputs: SL_insertion_sites.bed/.tsv, polyA_addition_sites.bed/.tsv,
    BAMs for SL and polyA, logs in outdir/logs/.
  - Delete a step marker (.stepX.done) or an output to re-run that step.
EOF
}

# --------------------------- Defaults -------------------------
READS="/home/block/Desktop/02-PTUTcruzi/data/Tcruzi_ystrain.fastq.gz"
GENOME="/home/block/Desktop/02-PTUTcruzi/final/ref_assembly_genome/yc6/yc6.fna"
OUTDIR="/home/block/Desktop/02-PTUTcruzi/final/02_mature_transcripts"
PTU_GFF=""
THREADS=12
MAP_PRESET="map-ont"

# 39-nt SL sequence (user-provided)
SL_SEQ="AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG"
SL_MIN_OVERLAP=20
SL_MAX_ERROR=0.22
POLYA_MIN_RUN=15

CUTADAPT_ENV="cutadapt"
MINIMAP2_ENV="minimap2"
SAMTOOLS_ENV="samtools"
BEDTOOLS_ENV="bedtools_env"
SEQKIT_ENV="seqkit"  # optional; not required

# --------------------------- CLI Parse ------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --reads) READS="$2"; shift 2 ;;
    --genome) GENOME="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --ptu-gff) PTU_GFF="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --sl-seq) SL_SEQ="$2"; shift 2 ;;
    --sl-min-overlap) SL_MIN_OVERLAP="$2"; shift 2 ;;
    --sl-max-error) SL_MAX_ERROR="$2"; shift 2 ;;
    --polya-min-run) POLYA_MIN_RUN="$2"; shift 2 ;;
    --cutadapt-env) CUTADAPT_ENV="$2"; shift 2 ;;
    --minimap2-env) MINIMAP2_ENV="$2"; shift 2 ;;
    --samtools-env) SAMTOOLS_ENV="$2"; shift 2 ;;
    --bedtools-env) BEDTOOLS_ENV="$2"; shift 2 ;;
    *) fail "Unknown argument: $1 (use --help)";;
  esac
done

# --------------------------- Prep -----------------------------
mkdir -p "${OUTDIR}"/{tmp,logs}
TMPDIR="${OUTDIR}/tmp"

# Conda bootstrap
_conda_base="$(conda info --base 2>/dev/null || true)"
[[ -z "${_conda_base}" || ! -d "${_conda_base}" ]] && fail "Conda base not found. Is conda installed?"
# shellcheck disable=SC1091
source "${_conda_base}/etc/profile.d/conda.sh"

activate() {
  local env_name="$1"
  msg "Activating conda env: ${env_name}"
  conda activate "${env_name}" || fail "conda activate ${env_name} failed"
}

# Helpers
run_if_missing() {
  local outfile="$1"; shift
  local env="$1"; shift
  if [[ -s "${outfile}" ]]; then
    msg "[SKIP] ${outfile} exists."
  else
    msg "[RUN ] Producing ${outfile}"
    activate "${env}"
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
    bash -euo pipefail -c "$(cat)" || fail "Block failed: ${desc}"
    : > "${marker}"
  fi
}

# ---------------------- Genome index (mmi) --------------------
GENOME_MMI="${GENOME}.mmi"
if [[ ! -s "${GENOME_MMI}" ]]; then
  activate "${MINIMAP2_ENV}"
  msg "Building minimap2 index: ${GENOME_MMI}"
  minimap2 -d "${GENOME_MMI}" "${GENOME}"
else
  msg "[SKIP] Genome index exists: ${GENOME_MMI}"
fi

# ---------------------- Build sl.fa (with RC) ----------------
SL_FA="${OUTDIR}/sl.fa"
if [[ ! -s "${SL_FA}" ]]; then
  msg "Creating SL fasta (with reverse-complement): ${SL_FA}"
  RC_SL=$(printf "%s\n" "${SL_SEQ}" | awk '
    BEGIN{
      map["A"]="T"; map["T"]="A"; map["C"]="G"; map["G"]="C";
      map["a"]="t"; map["t"]="a"; map["c"]="g"; map["g"]="c";
    }
    { s=$0; out="";
      for(i=length(s); i>=1; i--){
        b=substr(s,i,1); out = out ((b in map)? map[b] : b);
      } print out }')
  {
    echo ">SL";    echo "${SL_SEQ}"
    echo ">SL_rc"; echo "${RC_SL}"
  } > "${SL_FA}"
else
  msg "[SKIP] SL fasta exists: ${SL_FA}"
fi

# ------------------- Step 1: SL with cutadapt ----------------
SL_TRIM_FQ="${OUTDIR}/reads_with_SL.trimmed.fastq.gz"
SL_IDS="${OUTDIR}/reads_with_SL.ids"
run_block "${OUTDIR}/.step1.done" "${CUTADAPT_ENV}" "Detecting SL-containing reads (cutadapt)" <<EOF
cutadapt \
  -b "file:${SL_FA}" \
  --overlap ${SL_MIN_OVERLAP} \
  -e ${SL_MAX_ERROR} \
  --trimmed-only \
  -o "${SL_TRIM_FQ}" \
  "${READS}" \
  > "${OUTDIR}/logs/cutadapt_SL.log" 2>&1

# Extract read IDs (kept after trimming)
if [[ "${SL_TRIM_FQ}" == *.gz ]]; then
  zcat "${SL_TRIM_FQ}" | awk 'NR%4==1 {sub(/^@/,""); print \$1}' > "${SL_IDS}"
else
  awk 'NR%4==1 {sub(/^@/,""); print \$1}' "${SL_TRIM_FQ}" > "${SL_IDS}"
fi
EOF

[[ -s "${SL_TRIM_FQ}" ]] || msg "[WARN] No SL-positive reads detected. Check SL sequence/params."

# -------- Step 2: Map SL-trimmed (no pipes) + sort/index -----
SL_SAM="${OUTDIR}/sl.trimmed.toGenome.sam"
SL_BAM="${OUTDIR}/sl.trimmed.toGenome.bam"
SL_SORT_BAM="${OUTDIR}/sl.trimmed.toGenome.sorted.bam"
SL_SORT_BAI="${SL_SORT_BAM}.bai"

if [[ -s "${SL_TRIM_FQ}" ]]; then
  run_if_missing "${SL_SAM}" "${MINIMAP2_ENV}" \
    minimap2 -t "${THREADS}" -ax "${MAP_PRESET}" "${GENOME_MMI}" "${SL_TRIM_FQ}" -o "${SL_SAM}"
else
  msg "[WARN] ${SL_TRIM_FQ} empty. Creating empty SAM to keep idempotency."
  : > "${SL_SAM}"
fi

if [[ -s "${SL_SAM}" ]]; then
  head -n1 "${SL_SAM}" | grep -q '^@HD' || msg "[WARN] SAM header missing; file may be corrupted."
fi

if [[ -s "${SL_SAM}" ]]; then
  run_if_missing "${SL_BAM}" "${SAMTOOLS_ENV}" \
    samtools view -@ "${THREADS}" -b -F 4 -o "${SL_BAM}" "${SL_SAM}"
fi

if [[ -s "${SL_BAM}" ]]; then
  run_if_missing "${SL_SORT_BAM}" "${SAMTOOLS_ENV}" \
    samtools sort -@ "${THREADS}" -o "${SL_SORT_BAM}" "${SL_BAM}"
fi

if [[ -s "${SL_SORT_BAM}" ]]; then
  if [[ -s "${SL_SORT_BAI}" ]]; then
    msg "[SKIP] ${SL_SORT_BAI} exists."
  else
    activate "${SAMTOOLS_ENV}"
    samtools index "${SL_SORT_BAM}"
  fi
else
  msg "[WARN] Skipping index: ${SL_SORT_BAM} missing/empty."
fi

# -------- Step 3: Call SL insertion sites (BED + TSV) --------
SL_BED_RAW="${OUTDIR}/sl.alignments.bed"
SL_SITES_BED="${OUTDIR}/SL_insertion_sites.bed"
SL_SITES_TSV="${OUTDIR}/SL_insertion_sites.counts.tsv"

if [[ -s "${SL_SORT_BAM}" ]]; then
  run_if_missing "${SL_BED_RAW}" "${BEDTOOLS_ENV}" \
    bedtools bamtobed -i "${SL_SORT_BAM}" > "${SL_BED_RAW}"
else
  msg "[WARN] ${SL_SORT_BAM} empty; creating empty BED."
  : > "${SL_BED_RAW}"
fi

run_block "${OUTDIR}/.step3.done" "${BEDTOOLS_ENV}" "Computing SL insertion sites (5' ends)" <<'EOF'
# 1-bp site at the 5' aligned end (strand-aware)
awk 'BEGIN{OFS="\t"}{
  chrom=$1; start=$2; end=$3; name=$4; score=($5==""?0:$5); strand=$6;
  if (strand=="+") { s=start; e=start+1 }
  else if (strand=="-") { s=end-1; e=end }
  else { next }
  print chrom, s, e, name, score, strand
}' "${SL_BED_RAW}" > "${SL_SITES_BED}"

# Collapse identical sites and count
awk 'BEGIN{OFS="\t"}{key=$1 FS $2 FS $3 FS $6; c[key]++}END{
  print "chrom","start","end","count","strand";
  for(k in c){split(k,a,FS); print a[1],a[2],a[3],c[k],a[4]}
}' "${SL_SITES_BED}" \
| sort -k1,1 -k2,2n > "${SL_SITES_TSV}"
EOF

# ------------- Step 4: polyA detection + mapping --------------
POLYA_TRIM_FQ="${OUTDIR}/reads_with_polyA.trimmed.fastq.gz"
POLYA_SAM="${OUTDIR}/polyA.trimmed.toGenome.sam"
POLYA_BAM="${OUTDIR}/polyA.trimmed.toGenome.bam"
POLYA_SORT_BAM="${OUTDIR}/polyA.trimmed.toGenome.sorted.bam"
POLYA_SORT_BAI="${POLYA_SORT_BAM}.bai"

run_if_missing "${POLYA_TRIM_FQ}" "${CUTADAPT_ENV}" \
  cutadapt -a "A{${POLYA_MIN_RUN}}" --trimmed-only -o "${POLYA_TRIM_FQ}" "${READS}" \
  > "${OUTDIR}/logs/cutadapt_polyA.log" 2>&1

if [[ -s "${POLYA_TRIM_FQ}" ]]; then
  run_if_missing "${POLYA_SAM}" "${MINIMAP2_ENV}" \
    minimap2 -t "${THREADS}" -ax "${MAP_PRESET}" "${GENOME_MMI}" "${POLYA_TRIM_FQ}" -o "${POLYA_SAM}"
else
  msg "[WARN] ${POLYA_TRIM_FQ} empty. Creating empty SAM."
  : > "${POLYA_SAM}"
fi

if [[ -s "${POLYA_SAM}" ]]; then
  run_if_missing "${POLYA_BAM}" "${SAMTOOLS_ENV}" \
    samtools view -@ "${THREADS}" -b -F 4 -o "${POLYA_BAM}" "${POLYA_SAM}"
fi

if [[ -s "${POLYA_BAM}" ]]; then
  run_if_missing "${POLYA_SORT_BAM}" "${SAMTOOLS_ENV}" \
    samtools sort -@ "${THREADS}" -o "${POLYA_SORT_BAM}" "${POLYA_BAM}"
fi

if [[ -s "${POLYA_SORT_BAM}" ]]; then
  if [[ -s "${POLYA_SORT_BAI}" ]]; then
    msg "[SKIP] ${POLYA_SORT_BAI} exists."
  else
    activate "${SAMTOOLS_ENV}"
    samtools index "${POLYA_SORT_BAM}"
  fi
else
  msg "[WARN] Skipping index: ${POLYA_SORT_BAM} missing/empty."
fi

# --------- Step 5: polyA addition sites (BED + TSV) -----------
POLYA_BED_RAW="${OUTDIR}/polyA.alignments.bed"
POLYA_SITES_BED="${OUTDIR}/polyA_addition_sites.bed"
POLYA_SITES_TSV="${OUTDIR}/polyA_addition_sites.counts.tsv"

if [[ -s "${POLYA_SORT_BAM}" ]]; then
  run_if_missing "${POLYA_BED_RAW}" "${BEDTOOLS_ENV}" \
    bedtools bamtobed -i "${POLYA_SORT_BAM}" > "${POLYA_BED_RAW}"
else
  msg "[WARN] ${POLYA_SORT_BAM} empty; creating empty BED."
  : > "${POLYA_BED_RAW}"
fi

run_block "${OUTDIR}/.step5.done" "${BEDTOOLS_ENV}" "Computing polyA addition sites (3' ends)" <<'EOF'
awk 'BEGIN{OFS="\t"}{
  chrom=$1; start=$2; end=$3; name=$4; score=($5==""?0:$5); strand=$6;
  if (strand=="+") { s=end-1; e=end }
  else if (strand=="-") { s=start; e=start+1 }
  else { next }
  print chrom, s, e, name, score, strand
}' "${POLYA_BED_RAW}" > "${POLYA_SITES_BED}"

awk 'BEGIN{OFS="\t"}{key=$1 FS $2 FS $3 FS $6; c[key]++}END{
  print "chrom","start","end","count","strand";
  for(k in c){split(k,a,FS); print a[1],a[2],a[3],c[k],a[4]}
}' "${POLYA_SITES_BED}" \
| sort -k1,1 -k2,2n > "${POLYA_SITES_TSV}"
EOF

# -------------------------- Summaries ------------------------
activate "${SAMTOOLS_ENV}" || true
TOTAL_READS=$( ( [[ "${READS}" == *.gz ]] && zcat "${READS}" || cat "${READS}" ) | awk 'END{print NR/4}' )
SL_READS=$( ( [[ -s "${SL_TRIM_FQ}" && "${SL_TRIM_FQ}" == *.gz ]] && zcat "${SL_TRIM_FQ}" || cat "${SL_TRIM_FQ}" 2>/dev/null || true ) | awk 'END{print (NR?NR/4:0)}' )
POLYA_READS=$( ( [[ -s "${POLYA_TRIM_FQ}" && "${POLYA_TRIM_FQ}" == *.gz ]] && zcat "${POLYA_TRIM_FQ}" || cat "${POLYA_TRIM_FQ}" 2>/dev/null || true ) | awk 'END{print (NR?NR/4:0)}' )

SL_SITES=$(awk 'NR>1' "${SL_SITES_TSV}" 2>/dev/null | wc -l | awk '{print $1}')
POLYA_SITES=$(awk 'NR>1' "${POLYA_SITES_TSV}" 2>/dev/null | wc -l | awk '{print $1}')

cat <<SUM
=== Summary ===
Total reads:                 ${TOTAL_READS}
Reads with SL motif:         ${SL_READS}
Unique SL insertion sites:   ${SL_SITES}
Reads with polyA tail:       ${POLYA_READS}
Unique polyA addition sites: ${POLYA_SITES}
Logs:
  - ${OUTDIR}/logs/cutadapt_SL.log
  - ${OUTDIR}/logs/cutadapt_polyA.log
SUM

# --------------- Optional: SL vs PTU orientation --------------
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
  msg "[INFO] PTU_GFF not provided — skipping orientation check."
fi

msg "[DONE] Outputs in: ${OUTDIR}"
msg "  - ${SL_SITES_TSV} (SL insertion sites with counts)"
msg "  - ${POLYA_SITES_TSV} (polyA addition sites with counts)"
msg "  - ${SL_SITES_BED} , ${POLYA_SITES_BED} (1-bp, strand-aware BED)"
