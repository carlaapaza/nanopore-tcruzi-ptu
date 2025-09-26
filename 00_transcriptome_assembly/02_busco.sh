#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Script: busco.sh
# Purpose: Run BUSCO on a genome-guided transcriptome assembly
# Mode: transcriptome
# Author: Carla Apaza (T. cruzi project)
# ============================================================

# ----------------------- Paths (edit) -----------------------
FASTA="../assembly/transcripts.genome_guided.fa"
LINEAGE="../busco_downloads/lineages/trypanosoma_odb12"
OUT="../busco"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate busco

msg() { echo -e "[$(date +'%F %T')] $*"; }

# ------------------------- Run BUSCO ------------------------
msg "[INFO] Starting BUSCO in transcriptome mode..."

busco \
  -i "${FASTA}" \
  -l "${LINEAGE}" \
  -m transcriptome \
  -o "${OUT##*/}" \
  --out_path "$(dirname "${OUT}")" \
  --cpu 8 \
  --offline

conda deactivate busco

msg "[INFO] BUSCO finished successfully."
