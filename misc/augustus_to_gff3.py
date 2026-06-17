#!/usr/bin/env python
#
# 2020-09-03 previous version
# 2022-05-02 fix for old format without tags for transcript features
# 2023-01-26 remove Geneious formatting
# 2026-06-04 revamp

'''
  augustus_to_gff3.py  v2.0 last modified 2026-06-04
  generate true GFF3 from augustus output
  involves creating mRNA features from transcript
  and exon from CDS
  also will formalize some ID-Parent tags

augustus_to_gff3.py augustus.gff > augustus_no_comments.gff3

  some formats include exons already, rather than needing to create them from CDS

'''
#!/usr/bin/env python3

import sys
from collections import defaultdict

EXON_LIKE = {"initial", "internal", "terminal", "single"}
UTR_MAP = { "5'-UTR": "five_prime_UTR", "3'-UTR": "three_prime_UTR" }
PASS_FEATURES = {"start_codon", "stop_codon", "tss", "tts", "intron"}

example_metadata = """
# Evidence for and against this transcript:
# % of transcript supported by hints (any source): 91.5
# CDS exons: 19/21
#      E:  19 
# CDS introns: 18/20
#      E:  18 
# 5'UTR exons and introns: 5/5
#      E:   5 
# 3'UTR exons and introns: 1/1
#      E:   1 
# hint groups fully obeyed: 7
#      E:   7 (trinity_c25331_g6_i3)
# incompatible hint groups: 6
#      E:   6 (trinity_c18396_g1_i1,trinity_c18396_g1_i3,trinity_c26868_g1_i3,...)
# end gene g4
"""

def parse_attrs(attributes):
    attributes = attributes.replace("GFF Attribute=", "")

    if "codon_start" in attributes:
        attributes = attributes.split(";", 1)[1]

    d = {}
    for field in attributes.split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
        elif " " in field:
            k, v = field.split(" ", 1)
        else:
            continue
        d[k.strip()] = v.strip().replace('"', "")

    return d, attributes


def get_transcript_id(attrd, rawattrs):
    if attrd.get("transcript_id"):
        return attrd["transcript_id"]
    if attrd.get("Parent"):
        return attrd["Parent"]
    if attrd.get("ID"):
        return attrd["ID"]
    return rawattrs.strip()


def sanitize_attr_key(key):
    return (
        key.strip()
        .replace("%", "pct")
        .replace("'", "")
        .replace('"', "")
        .replace("(", "")
        .replace(")", "")
        .replace("/", "_of_")
        .replace("-", "_")
        .replace(" ", "_")
    )


def sanitize_attr_value(value):
    return (
        value.strip()
        .replace(";", ",")
        .replace("=", ":")
        .replace("\t", " ")
        .replace("\n", " ")
        .replace("(", "")
        .replace(")", "")
        .replace(",...", "")
    )


def parse_gene_comment_metadata(comment_lines, keep_protein=False):
    meta = {}
    protein_parts = []
    in_protein = False

    for line in comment_lines:
        text = line.lstrip("#").strip()

        if text.startswith("protein sequence ="):
            in_protein = True
            seq = text.split("=", 1)[1].strip()
            seq = seq.lstrip("[").rstrip("]")
            protein_parts.append(seq)
            continue

        if in_protein:
            if text.endswith("]"):
                protein_parts.append(text.rstrip("]"))
                in_protein = False
                continue
            elif text.startswith("Evidence for"):
                in_protein = False
            else:
                protein_parts.append(text)
                continue

        if ":" in text:
            key, value = text.split(":", 1)
            key = sanitize_attr_key(key)
            value = sanitize_attr_value(value)
            if key and value:
                meta[key] = value

    if keep_protein and protein_parts:
        meta["protein_sequence"] = sanitize_attr_value("".join(protein_parts))

    return meta


def merge_intervals(intervals):
    intervals = sorted(intervals)
    if not intervals:
        return []

    merged = [list(intervals[0])]

    for start, end in intervals[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])

    return [(s, e) for s, e in merged]


def flush_gene_block(geneid, records, comment_lines, keep_protein=False):
    extra_meta = parse_gene_comment_metadata(comment_lines, keep_protein=keep_protein)

    top_records = []
    subfeature_records = defaultdict(list)
    exon_parts = defaultdict(list)

    for lsplits in records:
        feature = lsplits[2]
        attrd, rawattrs = parse_attrs(lsplits[8])

        if feature == "gene":
            if rawattrs and rawattrs != ".":
                geneid = rawattrs

            attrs = [f"ID={geneid}", f"Name={geneid}"]
            attrs.extend(f"{k}={v}" for k, v in extra_meta.items())
            lsplits[8] = ";".join(attrs)
            top_records.append(lsplits)

        elif feature == "transcript":
            transcriptid = get_transcript_id(attrd, rawattrs)
            lsplits[2] = "mRNA"
            lsplits[8] = f"ID={transcriptid};Parent={geneid};Name={transcriptid}"
            top_records.append(lsplits)

        elif feature == "CDS":
            transcriptid = get_transcript_id(attrd, rawattrs)
            lsplits[8] = f"ID={transcriptid}.cds;Parent={transcriptid}"
            subfeature_records[transcriptid].append(lsplits)

        elif feature in EXON_LIKE:
            transcriptid = get_transcript_id(attrd, rawattrs)
            exon_parts[transcriptid].append((int(lsplits[3]), int(lsplits[4])))

        elif feature in UTR_MAP:
            transcriptid = get_transcript_id(attrd, rawattrs)

            lsplits[2] = UTR_MAP[feature]
            lsplits[7] = "."
            lsplits[8] = f"Parent={transcriptid}"
            subfeature_records[transcriptid].append(lsplits)

            exon_parts[transcriptid].append((int(lsplits[3]), int(lsplits[4])))

        elif feature in PASS_FEATURES:
            transcriptid = get_transcript_id(attrd, rawattrs)
            lsplits[8] = f"Parent={transcriptid}"
            subfeature_records[transcriptid].append(lsplits)

    for rec in top_records:
        print("\t".join(rec))

        if rec[2] != "mRNA":
            continue

        transcriptid = rec[8].split("ID=", 1)[1].split(";", 1)[0]

        for i, (start, end) in enumerate(merge_intervals(exon_parts[transcriptid]), 1):
            exon = list(rec)
            exon[2] = "exon"
            exon[3] = str(start)
            exon[4] = str(end)
            exon[7] = "."
            exon[8] = f"ID={transcriptid}.exon{i};Parent={transcriptid}"
            print("\t".join(exon))

        for subrec in subfeature_records[transcriptid]:
            print("\t".join(subrec))


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: augustus_to_gff3.py augustus.gff > augustus_no_comments.gff3")

    infile = sys.argv[1]
    keep_protein = "--keep-protein" in sys.argv[2:]

    current_gene = None
    block_records = []
    block_comment_lines = []

    with open(infile) as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith("# start gene"):
                current_gene = line.strip().split()[-1]
                block_records = []
                block_comment_lines = []
                continue

            if line.startswith("# end gene"):
                end_gene = line.strip().split()[-1]

                if current_gene is None:
                    current_gene = end_gene

                flush_gene_block(
                    current_gene,
                    block_records,
                    block_comment_lines,
                    keep_protein=keep_protein,
                )

                current_gene = None
                block_records = []
                block_comment_lines = []
                continue

            if line.startswith("#"):
                if current_gene is not None:
                    block_comment_lines.append(line)
                continue

            if not line:
                continue

            lsplits = line.split("\t")
            if len(lsplits) < 9:
                continue

            if current_gene is not None:
                block_records.append(lsplits)

    if current_gene is not None and block_records:
        flush_gene_block(
            current_gene,
            block_records,
            block_comment_lines,
            keep_protein=keep_protein,
        )


if __name__ == "__main__":
    main()
