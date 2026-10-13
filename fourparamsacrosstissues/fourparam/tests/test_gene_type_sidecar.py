"""The gene_type side-car is read as `text[gi]`, so its geometry is load-bearing.

`docs/genetype/*.txt` is one character per gene in genes.tsv order; a file one
gene short or out of step silently gives every later gene its neighbor's type.
The published file is only checked when it exists.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import build_gene_type as bg

REPO = Path(__file__).resolve().parents[3]
DOCS = REPO / "docs"
MANIFEST = DOCS / "manifest.json"

GTF = [
    "##description: fake\n",
    'chr1\tHAVANA\tgene\t1\t9\t.\t+\t.\tgene_id "ENSG001.3"; gene_type "protein_coding"; gene_name "A";\n',
    'chr1\tHAVANA\ttranscript\t1\t9\t.\t+\t.\tgene_id "ENSG001.3"; gene_type "protein_coding";\n',
    'chr1\tHAVANA\tgene\t1\t9\t.\t+\t.\tgene_id "ENSG002.1"; gene_type "lncRNA"; gene_name "B";\n',
    'chr1\tHAVANA\tgene\t1\t9\t.\t+\t.\tgene_id "ENSG003.7"; gene_type "IG_V_pseudogene";\n',
]


def test_gtf_reads_gene_records_only():
    lut = bg.read_gtf_gene_types(GTF)
    assert lut == {"ENSG001": "protein_coding", "ENSG002": "lncRNA",
                   "ENSG003": "IG_V_pseudogene"}


def test_encode_round_trips():
    order = ["ENSG002", "ENSG009", "ENSG001", "ENSG003"]
    text, classes = bg.encode(order, bg.read_gtf_gene_types(GTF))
    assert len(text) == len(order)
    assert text[1] == bg.MISSING
    assert [classes.get(c) for c in text] == \
        ["lncRNA", None, "protein_coding", "IG_V_pseudogene"]


def test_groups():
    assert bg.group_of("protein_coding") == "protein_coding"
    assert bg.group_of("lncRNA") == "lncRNA"
    assert bg.group_of("miRNA") == "small_ncRNA"
    # A pseudogene of an immune gene is a pseudogene, not an immune gene.
    assert bg.group_of("IG_V_pseudogene") == "pseudogene"
    assert bg.group_of("TR_V_gene") == "IG_TR"
    assert bg.group_of("TEC") == "other"


def test_published_file_matches_gene_index():
    if not MANIFEST.exists():
        return
    gt = json.loads(MANIFEST.read_text(encoding="utf-8")).get("gene_type")
    if not gt:
        return
    text = (DOCS / gt["base_url"] / gt["file"]).read_text(encoding="ascii")
    n_genes = sum(1 for line in (DOCS / "genes.tsv").open(encoding="utf-8")
                  if line.strip())
    assert len(text) == n_genes
    assert set(text) <= set(gt["classes"]) | {gt["missing_char"]}
    assert set(gt["groups"]) == set(gt["classes"].values())
