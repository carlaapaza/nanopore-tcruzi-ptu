#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1_annotate_ptus_ssr.py

Hybrid PTU/SSR annotator for trypanosomatid genomes.

Combines:
- Robust file handling & (optional) GFF3 header preservation like the "alexranieri" script
- Clean outputs (TSV / GTF / BED) like your original script
- Optional TSS emission at dSSRs (configurable)
- Works with GFF3 (key=value) or GTF (key "value") attributes

DEFINITIONS
- PTU: continuous run of selected "gene-like" features on the same strand within a contig
- dSSR (divergent): boundary/interval between a negative-strand PTU (←) and the next positive-strand PTU (→)
- cSSR (convergent): boundary/interval between a positive-strand PTU (→) and the next negative-strand PTU (←)

OUTPUTS (with -o PREFIX):
  PREFIX.ptu.tsv : chrom, ptu_id, start, end, strand, n_genes, gene_ids
  PREFIX.ptu.gtf : 'ptu' features with attributes ID, Name, n_genes, gene_ids
  PREFIX.ptu.bed : BED6 for PTUs

  PREFIX.ssr.tsv : chrom, ssr_id, start, end, type, left_ptu_id, right_ptu_id
  PREFIX.ssr.gtf : 'SSR' features with attributes ID, Type, LeftPTU, RightPTU
  PREFIX.ssr.bed : BED6 for SSRs (score=0, strand='.')

OPTIONAL:
  PREFIX.annotated.gff3 : original headers + appended PTU/SSR/TSS features (GFF3 style)
  TSS features at dSSRs (window centered on the SSR; width = fraction * SSR length)

USAGE
  python3 1_annotate_ptus_ssr.py \
    -i /home/block/Desktop/02-PTUTcruzi/final/1_references/gff/TriTrypDB-68_TcruziYC6.gff \
    -o /home/block/Desktop/02-PTUTcruzi/final/2_ptu_ssr_annotation/yc6 \
    --features CDS \
    --emit-tss --tss-fraction 0.5 \
    --append-gff3
"""

import argparse
import sys
from collections import defaultdict

def parse_args():
    ap = argparse.ArgumentParser(description="Detect PTUs and SSRs from a GFF3/GTF.")
    ap.add_argument("-i", "--gff", required=True, help="Input GFF3/GTF file.")
    ap.add_argument("-o", "--out-prefix", required=True, help="Output prefix for TSV/GTF/BED files.")
    # Feature(s) that define “genes” for PTU building (comma-separated)
    ap.add_argument("--features", default="gene", 
                    help="Comma-separated feature types to treat as genes (e.g. 'transcript' or 'gene,transcript').")
    ap.add_argument("--include-biotype", action="append", default=None,
                    help="If set, only keep entries whose attribute gene_biotype is in this list. "
                         "Use multiple --include-biotype flags to add more.")
    ap.add_argument("--append-gff3", action="store_true",
                    help="Write PREFIX.annotated.gff3 with original headers + appended PTU/SSR/TSS features.")
    ap.add_argument("--emit-tss", action="store_true",
                    help="Emit TSS windows at dSSRs (like alexranieri script).")
    ap.add_argument("--tss-fraction", type=float, default=0.5,
                    help="TSS window width as a fraction of SSR length (default 0.5; the window is centered).")
    return ap.parse_args()

# ---------- Attribute parsing: supports GFF3 and GTF ----------

def parse_attrs(attr_field):
    """
    Parse GFF3 (key=value;key2=value2) or GTF (key "value"; key2 "value2";) attributes into a dict.
    Handles empty attributes gracefully.
    """
    d = {}
    if not attr_field or attr_field == ".":
        return d

    s = attr_field.strip()

    # Heuristic: if we see '=' it's likely GFF3; otherwise treat as GTF-like
    if "=" in s and ";" in s and not ('"' in s and s.count('"') % 2 == 0 and s.find("=") > s.find('"')):
        # GFF3 style: key=value;key2=value2
        parts = [p.strip() for p in s.split(";") if p.strip()]
        for p in parts:
            if "=" in p:
                key, val = p.split("=", 1)
                d[key.strip()] = val.strip()
            else:
                # tolerate bare tokens
                d[p.strip()] = ""
    else:
        # GTF style: key "value"; key2 "value2";
        parts = [p.strip() for p in s.split(";") if p.strip()]
        for p in parts:
            if " " in p:
                key, val = p.split(" ", 1)
                d[key.strip()] = val.strip().strip('"')
            else:
                d[p.strip()] = ""

    return d

# ---------- I/O and record handling ----------

def read_records(gff_path):
    """
    Read all non-comment records from GFF/GTF. Keep header lines separately if needed.
    Return (headers:list[str], records:list[dict])
    Each record: dict with keys chrom, source, feature, start, end, score, strand, frame, attrs(dict), raw(list[str])
    """
    headers = []
    records = []
    with open(gff_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                headers.append(line.rstrip("\n"))
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = cols
            try:
                start_i = int(start); end_i = int(end)
            except ValueError:
                continue
            # Normalize strand
            if strand not in ("+","-","."):
                strand = "."
            rec = {
                "chrom": chrom, "source": source, "feature": feature,
                "start": start_i, "end": end_i, "score": score,
                "strand": strand, "frame": frame,
                "attrs": parse_attrs(attrs),
                "raw": cols
            }
            records.append(rec)
    # sort by chrom, start, end
    records.sort(key=lambda r: (r["chrom"], r["start"], r["end"], r["feature"]))
    return headers, records

def filter_gene_like(records, feature_set, include_biotypes=None):
    """
    Keep only records whose feature is in feature_set and (if provided) gene_biotype is allowed.
    Return dict: chrom -> list of gene-like dicts (chrom,start,end,strand,gene_id,attrs)
    """
    genes_by_chrom = defaultdict(list)
    allowed = set(include_biotypes) if include_biotypes else None

    for r in records:
        if r["feature"] not in feature_set:
            continue
        ad = r["attrs"]
        if allowed is not None:
            gb = ad.get("gene_biotype") or ad.get("biotype") or ad.get("transcript_biotype")
            if gb not in allowed:
                continue
        # Pick a stable gene ID
        gid = (ad.get("gene_id") or ad.get("locus_tag") or ad.get("ID") or
               ad.get("Name") or ad.get("orig_transcript_id") or
               ad.get("orig_protein_id") or f"{r['chrom']}:{r['start']}-{r['end']}:{r['strand']}")
        if r["strand"] not in ("+","-"):
            # skip unstranded for PTU logic
            continue
        genes_by_chrom[r["chrom"]].append({
            "chrom": r["chrom"],
            "start": r["start"],
            "end": r["end"],
            "strand": r["strand"],
            "gene_id": gid,
            "attrs": ad
        })

    # sort each chrom by start
    for c in genes_by_chrom:
        genes_by_chrom[c].sort(key=lambda g: (g["start"], g["end"]))
    return genes_by_chrom

# ---------- PTU / SSR core logic ----------

def build_ptus(genes_by_chrom):
    """
    From sorted genes, build PTUs as runs of same-strand genes.
    Returns list of PTU dicts with keys:
      chrom, start, end, strand, genes (list of gene_ids), n_genes, ptu_id, index
    PTU IDs are PTU:<chrom>:<1-based index per chrom>
    """
    ptus = []
    for chrom, genes in genes_by_chrom.items():
        if not genes:
            continue
        current = None
        idx = 0
        for g in genes:
            if current is None:
                idx += 1
                current = {
                    "chrom": chrom,
                    "start": g["start"],
                    "end": g["end"],
                    "strand": g["strand"],
                    "genes": [g["gene_id"]],
                    "n_genes": 1,
                    "index": idx,
                }
            else:
                if g["strand"] == current["strand"]:
                    current["end"] = max(current["end"], g["end"])
                    current["genes"].append(g["gene_id"])
                    current["n_genes"] += 1
                else:
                    current["ptu_id"] = f"PTU:{chrom}:{current['index']}"
                    ptus.append(current)
                    idx += 1
                    current = {
                        "chrom": chrom,
                        "start": g["start"],
                        "end": g["end"],
                        "strand": g["strand"],
                        "genes": [g["gene_id"]],
                        "n_genes": 1,
                        "index": idx,
                    }
        if current is not None:
            current["ptu_id"] = f"PTU:{chrom}:{current['index']}"
            ptus.append(current)
    return ptus

def derive_ssrs(ptus):
    """
    Create SSRs between adjacent PTUs on the same chromosome.
    Interval:
      start = prev.end + 1
      end   = next.start - 1
    If start > end => collapse to 1-bp boundary at min(prev_end, next_start)
    Type:
      dSSR if prev strand '-' and next strand '+'
      cSSR if prev strand '+' and next strand '-'
      'inter-PTU' otherwise (should be rare if PTUs built as runs)
    """
    by_chrom = defaultdict(list)
    for p in ptus:
        by_chrom[p["chrom"]].append(p)
    for c in by_chrom:
        by_chrom[c].sort(key=lambda x: (x["start"], x["end"]))
    ssrs = []
    d_count = 0
    c_count = 0
    inter_count = 0
    for chrom, plist in by_chrom.items():
        for i in range(len(plist) - 1):
            left = plist[i]
            right = plist[i+1]
            start = left["end"] + 1
            end = right["start"] - 1
            if start > end:
                boundary = min(left["end"], right["start"])
                start = boundary
                end = boundary
            lt = left["strand"]; rt = right["strand"]
            if lt == "-" and rt == "+":
                ssr_type = "dSSR"; d_count += 1
                ssr_id = f"dSSR_{d_count}"
            elif lt == "+" and rt == "-":
                ssr_type = "cSSR"; c_count += 1
                ssr_id = f"cSSR_{c_count}"
            else:
                ssr_type = "inter-PTU"; inter_count += 1
                ssr_id = f"IPR_{inter_count}"
            ssrs.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "type": ssr_type,
                "left_ptu_id": left["ptu_id"],
                "right_ptu_id": right["ptu_id"],
                "ssr_id": ssr_id
            })
    return ssrs

# ---------- Writers ----------

def gtf_attr_str(d):
    return " ".join(f'{k} "{v}";' for k, v in d.items())

def write_ptu_tsv(ptus, path):
    with open(path, "w", encoding="utf-8") as out:
        out.write("\t".join(["chrom","ptu_id","start","end","strand","n_genes","gene_ids"]) + "\n")
        for p in ptus:
            out.write("\t".join([
                p["chrom"], p["ptu_id"], str(p["start"]), str(p["end"]), p["strand"],
                str(p["n_genes"]), ",".join(p["genes"])
            ]) + "\n")

def write_ptu_bed(ptus, path):
    with open(path, "w", encoding="utf-8") as out:
        for p in ptus:
            bed_start = max(0, p["start"] - 1)  # BED is 0-based, half-open
            out.write("\t".join([
                p["chrom"], str(bed_start), str(p["end"]),
                p["ptu_id"], str(p["n_genes"]), p["strand"]
            ]) + "\n")

def write_ptu_gtf(ptus, path, source="PTUAnnotator"):
    with open(path, "w", encoding="utf-8") as out:
        out.write("##gff-version 2.2\n")
        for p in ptus:
            attrs = {
                "ID": p["ptu_id"],
                "Name": p["ptu_id"],
                "n_genes": str(p["n_genes"]),
                "gene_ids": ",".join(p["genes"])
            }
            row = [
                p["chrom"], source, "ptu",
                str(p["start"]), str(p["end"]),
                ".", p["strand"], ".", gtf_attr_str(attrs)
            ]
            out.write("\t".join(row) + "\n")

def write_ssr_tsv(ssrs, path):
    with open(path, "w", encoding="utf-8") as out:
        out.write("\t".join(["chrom","ssr_id","start","end","type","left_ptu_id","right_ptu_id"]) + "\n")
        for s in ssrs:
            out.write("\t".join([
                s["chrom"], s["ssr_id"], str(s["start"]), str(s["end"]),
                s["type"], s["left_ptu_id"], s["right_ptu_id"]
            ]) + "\n")

def write_ssr_bed(ssrs, path):
    with open(path, "w", encoding="utf-8") as out:
        for s in ssrs:
            bed_start = max(0, s["start"] - 1)
            out.write("\t".join([
                s["chrom"], str(bed_start), str(s["end"]),
                s["ssr_id"], "0", "."
            ]) + "\n")

def write_ssr_gtf(ssrs, path, source="PTUAnnotator"):
    with open(path, "w", encoding="utf-8") as out:
        out.write("##gff-version 2.2\n")
        for s in ssrs:
            attrs = {
                "ID": s["ssr_id"],
                "Type": s["type"],
                "LeftPTU": s["left_ptu_id"],
                "RightPTU": s["right_ptu_id"]
            }
            row = [
                s["chrom"], source, "SSR",
                str(s["start"]), str(s["end"]),
                ".", ".", ".", gtf_attr_str(attrs)
            ]
            out.write("\t".join(row) + "\n")

def write_annotated_gff3(headers, ptus, ssrs, path, emit_tss=False, tss_fraction=0.5):
    """
    Append our features in GFF3 style (key=value) after original headers.
    Adds TSS for dSSRs if emit_tss=True. TSS window centered with size=fraction*SSR_length.
    """
    with open(path, "w", encoding="utf-8") as out:
        for h in headers:
            out.write(h + "\n")
        prog = "annotatePolycistron"
        # PTUs
        cds = trna = ncrna = rrna = snorna = 0  # kept for compatibility naming if needed
        for p in ptus:
            # ID=PTU_<chrom>_<index>
            pid = p["ptu_id"].replace("PTU:", "PTU_").replace(":", "_")
            desc = f"ID={pid};contentCount={p['n_genes']};content={','.join(p['genes'])}"
            row = [p["chrom"], prog, "ptu",
                   str(p["start"]), str(p["end"]),
                   ".", p["strand"], ".", desc]
            out.write("\t".join(row) + "\n")
        # SSRs + TSS
        tss_counter = 0
        for s in ssrs:
            lpid = s["left_ptu_id"].replace("PTU:", "PTU_").replace(":", "_")
            rpid = s["right_ptu_id"].replace("PTU:", "PTU_").replace(":", "_")
            desc = f"ID={s['ssr_id']};adjacentPol={lpid},{rpid}"
            row = [s["chrom"], prog, s["type"], str(s["start"]), str(s["end"]), ".",
                   ".", ".", desc]
            out.write("\t".join(row) + "\n")

            if emit_tss and s["type"] == "dSSR":
                # build TSS window centered at SSR midpoint; size = fraction * SSR length (>=1)
                L = max(1, s["end"] - s["start"] + 1)
                size = max(1, int(round(L * float(tss_fraction))))
                mid = (s["start"] + s["end"]) // 2
                tss_start = max(1, mid - size // 2)
                tss_end = tss_start + size - 1
                tss_counter += 1
                tss_id = f"ID=TSS_{tss_counter};adjacentPol={lpid},{rpid}"
                tss_row = [s["chrom"], prog, "TSS",
                           str(tss_start), str(tss_end),
                           ".", ".", ".", tss_id]
                out.write("\t".join(tss_row) + "\n")

# ---------- MAIN ----------

def main():
    args = parse_args()

    feature_set = {f.strip() for f in args.features.split(",") if f.strip()}
    headers, records = read_records(args.gff)

    # select gene-like features
    genes_by_chrom = filter_gene_like(records, feature_set, include_biotypes=args.include_biotype)
    if not genes_by_chrom:
        sys.stderr.write("No gene-like features found. "
                         "Check --features and/or --include-biotype filters.\n")
        sys.exit(1)

    # build PTUs and SSRs
    ptus = build_ptus(genes_by_chrom)
    ptus.sort(key=lambda p: (p["chrom"], p["start"], p["end"]))
    ssrs = derive_ssrs(ptus)

    # write outputs
    write_ptu_tsv(ptus, f"{args.out_prefix}.ptu.tsv")
    write_ptu_bed(ptus, f"{args.out_prefix}.ptu.bed")
    write_ptu_gtf(ptus, f"{args.out_prefix}.ptu.gtf")

    write_ssr_tsv(ssrs, f"{args.out_prefix}.ssr.tsv")
    write_ssr_bed(ssrs, f"{args.out_prefix}.ssr.bed")
    write_ssr_gtf(ssrs, f"{args.out_prefix}.ssr.gtf")

    if args.append_gff3:
        write_annotated_gff3(headers, ptus, ssrs, f"{args.out_prefix}.annotated.gff3",
                             emit_tss=args.emit_tss, tss_fraction=args.tss_fraction)

    # Summary
    n_ptus = len(ptus)
    n_ssrs = len(ssrs)
    n_d = sum(1 for s in ssrs if s["type"] == "dSSR")
    n_c = sum(1 for s in ssrs if s["type"] == "cSSR")
    sys.stderr.write(f"PTUs: {n_ptus} | SSRs: {n_ssrs} (dSSR={n_d}, cSSR={n_c})\n")
    sys.stderr.write(f"Wrote: {args.out_prefix}.ptu.tsv/.gtf/.bed and {args.out_prefix}.ssr.tsv/.gtf/.bed\n")
    if args.append_gff3:
        sys.stderr.write(f"Wrote: {args.out_prefix}.annotated.gff3\n")

if __name__ == "__main__":
    main()
