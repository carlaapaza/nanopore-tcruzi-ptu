#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2_classify_ptus.py

Classify genes into genomic compartments (core, disruptive, both, other)
using annotation text prioritized from 'description' (GFF3) and 'product'
(compat), then classify PTUs with the 80% rule:

- "PTU core"        if >=80% of member genes are 'core'
- "PTU disruptive"  if >=80% are 'disruptive'
- "PTU both"        if >=80% are 'both' (GP63/RHS/DGF-1 bucket)
- otherwise "PTU mix"

Inputs:
  --genes-gtf : GFF3/GTF with gene-level entries (feature given by --gene-feature)
  --ptu-gtf   : GTF with PTU features (feature 'ptu')
  --out       : output prefix

Outputs:
  <out>.gene_comp.tsv : gene_id, chrom, start, end, strand, ptu_id, description, gene_compartment
  <out>.ptu_comp.tsv  : per-PTU counts, fractions, and PTU class
"""

import argparse
import csv
import re
import sys
from collections import defaultdict

# -------------------- args --------------------

def parse_args():
    ap = argparse.ArgumentParser(description="Classify genes and PTUs into genomic compartments and PTU classes.")
    ap.add_argument("--genes-gtf", required=True,
                    help="GFF3/GTF with gene entries (e.g., 'protein_coding_gene', 'transcript', 'gene', 'CDS').")
    ap.add_argument("--ptu-gtf", required=True,
                    help="PTU GTF with features named 'ptu'.")
    ap.add_argument("--gene-feature", default="protein_coding_gene",
                    help="Feature to treat as genes (default: protein_coding_gene).")
    ap.add_argument("--out", required=True, help="Output prefix.")
    ap.add_argument("--strand-match", action="store_true",
                    help="Require gene and PTU to be on same strand when assigning membership.")

    # Patterns (case-insensitive)
    ap.add_argument(
        "--disruptive-patterns",
        default=r"(\btrans[- ]?sialidase\b|\bts\d+\b|\bmucin\b|\bmasp\b|\bmucin[- ]?associated\b|\btc[- ]?muc\b)",
        help="Regex for disruptive multigene families (TS, mucin, MASP, TcMUC, mucin-associated)."
    )
    ap.add_argument(
        "--both-patterns",
        default=r"(\bgp ?63\b|dispersed gene family 1|retrotransposon hot spot|\brhs\b|\bdgf[- ]?1\b|\bdgf1\b)",
        help="Regex for 'both' bucket (GP63, Dispersed Gene Family 1, RHS, DGF-1)."
    )
    ap.add_argument(
        "--core-patterns",
        default=r"(hypothetical conserved|conserved|hypothetical protein)",
        help="Regex for core-like descriptors (includes 'hypothetical protein')."
    )

    # Prioridad de atributos para el texto de clasificación y la columna description
    # Para GFF3 tuyo: description e ebi_biotype son comunes; mantenemos product/note por compatibilidad.
    ap.add_argument(
        "--attr-fields",
        default="description,product,note,Name,gene_name,ebi_biotype,gene_biotype,locus_tag",
        help="Comma-separated attribute keys to search for patterns (priority)."
    )

    ap.add_argument("--ptu-pure-threshold", type=float, default=0.80,
                    help="Fraction threshold for PTU compartment call (default: 0.80).")

    # Permitir excluir 'hypothetical protein' del core si se desea (por defecto, incluido)
    ap.add_argument("--exclude-plain-hypothetical-from-core", action="store_true",
                    help="If set, 'hypothetical protein' WITHOUT 'conserved' won't count as core.")
    return ap.parse_args()

# -------------------- parsing helpers --------------------

def parse_attrs(attr_field):
    """
    Parse GFF3 (key=value;key2=value2) or GTF (key "value"; key2 "value2";) into dict.
    """
    d = {}
    if not attr_field or attr_field == ".":
        return d
    s = attr_field.strip()

    # Heurística: si hay '=' y ';', es más probable GFF3
    if "=" in s and ";" in s and not ('"' in s and s.count('"') % 2 == 0 and s.find("=") > s.find('"')):
        parts = [p.strip() for p in s.split(";") if p.strip()]
        for p in parts:
            if "=" in p:
                k, v = p.split("=", 1)
                d[k.strip()] = v.strip()
            else:
                d[p.strip()] = ""
    else:
        parts = [p.strip() for p in s.split(";") if p.strip()]
        for p in parts:
            if " " in p:
                k, v = p.split(" ", 1)
                d[k.strip()] = v.strip().strip('"')
            else:
                d[p.strip()] = ""
    return d

def read_gtf_like(path, feature_filter=None):
    feats = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            chrom, src, feat, start, end, score, strand, frame, attrs = f
            if feature_filter and feat != feature_filter:
                continue
            try:
                start_i = int(start); end_i = int(end)
            except ValueError:
                continue
            if strand not in ("+","-","."):
                strand = "."
            feats.append({
                "chrom": chrom,
                "source": src,
                "feature": feat,
                "start": start_i,
                "end": end_i,
                "strand": strand,
                "attrs": parse_attrs(attrs)
            })
    return feats

def get_first(d, keys, fallback=None):
    for k in keys:
        if k in d and d[k]:
            return d[k]
    return fallback

def overlap(a1, a2, b1, b2):
    return not (a2 < b1 or b2 < a1)

def build_index_by_chrom(feats):
    byc = defaultdict(list)
    for x in feats:
        byc[x["chrom"]].append(x)
    for c in byc:
        byc[c].sort(key=lambda r: (r["start"], r["end"]))
    return byc

# -------------------- classification --------------------

def classify_gene_text(text, re_disr, re_both, re_core, exclude_plain_hyp=False):
    """
    Return: 'disruptive' | 'both' | 'core' | 'other'
    """
    s = (text or "").lower()
    if not s:
        return "other"
    if re_disr.search(s):
        return "disruptive"
    if re_both.search(s):
        return "both"
    if exclude_plain_hyp and ("hypothetical protein" in s) and ("conserved" not in s):
        return "other"
    if re_core.search(s):
        return "core"
    return "other"

# -------------------- main --------------------

def main():
    args = parse_args()

    # Compilar patrones y campos
    re_disr = re.compile(args.disruptive_patterns, re.IGNORECASE)
    re_both = re.compile(args.both_patterns, re.IGNORECASE)
    re_core = re.compile(args.core_patterns, re.IGNORECASE)
    fields = [f.strip() for f in args.attr_fields.split(",") if f.strip()]

    # Cargar datos
    genes = read_gtf_like(args.genes_gtf, feature_filter=args.gene_feature)
    ptus  = read_gtf_like(args.ptu_gtf,   feature_filter="ptu")

    if not genes:
        sys.exit(f"No features found for --gene-feature '{args.gene_feature}' in {args.genes_gtf}.")
    if not ptus:
        sys.exit("No PTU features found (feature must be 'ptu').")

    # Normalizar PTU IDs
    for p in ptus:
        p["ptu_id"] = get_first(p["attrs"], ["ID","Name"], f"PTU:{p['chrom']}:{p['start']}-{p['end']}")

    # Indexar por cromosoma
    genes_byc = build_index_by_chrom(genes)
    ptus_byc  = build_index_by_chrom(ptus)

    # Asignar genes a PTUs + clasificar
    gene_rows = []
    classes_by_ptu = defaultdict(list)

    for chrom, plist in ptus_byc.items():
        if chrom not in genes_byc:
            continue
        glist = genes_byc[chrom]
        for p in plist:
            for g in glist:
                if args.strand_match and g["strand"] != p["strand"]:
                    continue
                if overlap(p["start"], p["end"], g["start"], g["end"]):
                    # construir texto para clasificación (description/product primero)
                    vals = []
                    for k in fields:
                        v = g["attrs"].get(k)
                        if v:
                            vals.append(v)
                    text = " | ".join(vals)
                    gene_class = classify_gene_text(
                        text, re_disr, re_both, re_core,
                        exclude_plain_hyp=args.exclude_plain_hypothetical_from_core
                    )

                    gene_id = get_first(
                        g["attrs"],
                        ["gene_id","ID","locus_tag","Name","orig_protein_id","orig_transcript_id"],
                        f"{chrom}:{g['start']}-{g['end']}"
                    )
                    description = (g["attrs"].get("description") or
                                   g["attrs"].get("product") or
                                   g["attrs"].get("note") or
                                   get_first(g["attrs"], ["Name","gene_name"], ""))

                    gene_rows.append({
                        "gene_id": gene_id,
                        "chrom": chrom,
                        "start": g["start"],
                        "end": g["end"],
                        "strand": g["strand"],
                        "ptu_id": p["ptu_id"],
                        "description": description,
                        "gene_compartment": gene_class
                    })
                    classes_by_ptu[p["ptu_id"]].append(gene_class)

    # Escribir gene_comp.tsv
    gene_tsv = f"{args.out}.gene_comp.tsv"
    with open(gene_tsv, "w", encoding="utf-8", newline="") as out:
        cols = ["gene_id","chrom","start","end","strand","ptu_id","description","gene_compartment"]
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in gene_rows:
            w.writerow(r)

    # Resumen por PTU + regla 80%
    ptu_rows = []
    ptu_lookup = {p["ptu_id"]: p for p in ptus}
    thr = args.ptu_pure_threshold

    for ptu_id, cls_list in classes_by_ptu.items():
        n = len(cls_list)
        n_core = sum(1 for c in cls_list if c == "core")
        n_disr = sum(1 for c in cls_list if c == "disruptive")
        n_both = sum(1 for c in cls_list if c == "both")
        n_other = n - n_core - n_disr - n_both

        f_core = (n_core / n) if n else 0.0
        f_disr = (n_disr / n) if n else 0.0
        f_both = (n_both / n) if n else 0.0

        if f_core >= thr:
            pclass = "PTU core"
        elif f_disr >= thr:
            pclass = "PTU disruptive"
        elif f_both >= thr:
            pclass = "PTU both"
        else:
            pclass = "PTU mix"

        p = ptu_lookup[ptu_id]
        ptu_rows.append({
            "ptu_id": ptu_id,
            "chrom": p["chrom"],
            "start": p["start"],
            "end": p["end"],
            "strand": p["strand"],
            "n_genes": n,
            "n_core": n_core,
            "n_disruptive": n_disr,
            "n_both": n_both,
            "n_other": n_other,
            "frac_core": f"{f_core:.3f}",
            "frac_disruptive": f"{f_disr:.3f}",
            "frac_both": f"{f_both:.3f}",
            "ptu_class": pclass
        })

    # PTUs sin genes asignados -> PTU mix
    covered = set(classes_by_ptu.keys())
    for p in ptus:
        if p["ptu_id"] in covered:
            continue
        ptu_rows.append({
            "ptu_id": p["ptu_id"],
            "chrom": p["chrom"],
            "start": p["start"],
            "end": p["end"],
            "strand": p["strand"],
            "n_genes": 0,
            "n_core": 0, "n_disruptive": 0, "n_both": 0, "n_other": 0,
            "frac_core": "0.000", "frac_disruptive": "0.000", "frac_both": "0.000",
            "ptu_class": "PTU mix"
        })

    ptu_rows.sort(key=lambda r: (r["chrom"], int(r["start"]), int(r["end"])))

    ptu_tsv = f"{args.out}.ptu_comp.tsv"
    with open(ptu_tsv, "w", encoding="utf-8", newline="") as out:
        cols = ["ptu_id","chrom","start","end","strand","n_genes",
                "n_core","n_disruptive","n_both","n_other",
                "frac_core","frac_disruptive","frac_both","ptu_class"]
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in ptu_rows:
            w.writerow(r)

    sys.stderr.write(f"[done] Wrote {gene_tsv}\n")
    sys.stderr.write(f"[done] Wrote {ptu_tsv}\n")

if __name__ == "__main__":
    main()
