#!/usr/bin/env python3
"""Render paper/paper.md with the embedded Times New Roman font from the user-owned font package.

The container's pdflatex is too minimal (missing geometry.sty etc.), so this script
renders directly. Figures from figures/ are appended as a captioned section.
"""
import os
import re

from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, Image, PageBreak)

FONTDIR = os.environ.get("MEGA27_TNR_DIR")
if not FONTDIR:
    raise RuntimeError("Set MEGA27_TNR_DIR to the folder containing user-owned Times New Roman TTFs")
for _name, _file in [("TNR", "Times.TTF"), ("TNR-Bold", "Timesbd.TTF"),
                     ("TNR-Italic", "Timesi.TTF"), ("TNR-BoldItalic", "Timesbi.TTF")]:
    pdfmetrics.registerFont(TTFont(_name, os.path.join(FONTDIR, _file)))
pdfmetrics.registerFontFamily("TNR", normal="TNR", bold="TNR-Bold",
                              italic="TNR-Italic", boldItalic="TNR-BoldItalic")
ROOT = os.path.join(os.path.dirname(__file__), "..")
STYLES = {
    "title": ParagraphStyle("title", fontName="TNR-Bold", fontSize=17, leading=22, spaceAfter=6),
    "author": ParagraphStyle("author", fontName="TNR", fontSize=12, leading=15, spaceAfter=2),
    "h1": ParagraphStyle("h1", fontName="TNR-Bold", fontSize=14, leading=18, spaceBefore=14, spaceAfter=6),
    "h2": ParagraphStyle("h2", fontName="TNR-Bold", fontSize=12, leading=15, spaceBefore=10, spaceAfter=4),
    "body": ParagraphStyle("body", fontName="TNR", fontSize=11.5, leading=16, spaceAfter=7, alignment=4),
    "caption": ParagraphStyle("caption", fontName="TNR-Italic", fontSize=10, leading=13, spaceBefore=4, spaceAfter=14),
    "cell": ParagraphStyle("cell", fontName="TNR", fontSize=8.5, leading=11),
    "cellb": ParagraphStyle("cellb", fontName="TNR-Bold", fontSize=8.5, leading=11),
}

FIGURES = [
    ("fig1_efficacy_distribution.png", "Figure 1. CNN efficacy scores of all NGG sites in the 110 kb dsx locus; red line marks the published Kyrou 2018 gRNA."),
    ("fig2_top_designs.png", "Figure 2. Top-20 v1 candidate sites (blue) and the published Kyrou gRNA at its true rank (red)."),
    ("fig3_drive_trajectories.png", "Figure 3. Drive allele frequency trajectories across homing rates (e = 0.01)."),
    ("fig4_resistance.png", "Figure 4. Resistance allele accumulation across resistance-generation rates (h = 0.99)."),
    ("fig5_wright_fisher.png", "Figure 5. Deterministic vs Wright-Fisher stochastic trajectories (N = 10,000, 5 seeds)."),
    ("fig6_constraint_landscape.png", "Figure 6. Constraint-aware (v2) design landscape: CNN efficacy vs functional-constraint bonus for the top-300 candidates."),
    ("fig7_genomewide_offtargets.png", "Figure 7. Exact genome-wide off-target burden of the v2 top-30 candidates and the published Kyrou guide, stacked by mismatch class (symlog scale; AgamP5, 246 Mb, <=3 mismatches)."),
    ("fig8_v3_ranking.png", "Figure 8. v3 composite ranking of the 28 filter-passing candidates against the published Kyrou 2018 guide (dashed line, v3 = 0.837). Red: named lead candidates dsx-v3-1 and dsx-v3-2."),
    ("fig9_locus_map.png", "Figure 9. Design map of the dsx locus: v3 candidate protospacers (teal pass, grey fail, red named leads), the published Kyrou 2018 guide (magenta, below axis) and the female-specific exon (dashed)."),
    ("fig10_enumerator_benchmark.png", "Figure 10. Enumerator throughput vs published Cas-OFFinder numbers, normalized to seconds per guide-Gb (log scale). Our enumerator is 523x slower than Cas-OFFinder CPU and 10,431x slower than GPU; see Appendix M for the honest reading."),
]


def build() -> str:
    src = open(os.path.join(ROOT, "paper", "paper.md")).read()
    src = re.sub(r"^---.*?^---\s*", "", src, flags=re.S | re.M)
    lines = src.split("\n")
    story, i = [], 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("# "):
            story.append(Paragraph(line[2:], STYLES["h1"]))
        elif line.startswith("## "):
            story.append(Paragraph(line[3:], STYLES["h2"]))
        elif line.startswith("|") and i + 1 < len(lines) and set(lines[i + 1]) <= set("|-: "):
            header = [c.strip() for c in line.strip("|").split("|")]
            rows = []
            i += 2
            while i < len(lines) and lines[i].startswith("|"):
                rows.append([c.strip() for c in lines[i].strip("|").split("|")])
                i += 1
            i -= 1
            data = [[Paragraph(c, STYLES["cellb"]) for c in header]] + \
                   [[Paragraph(c, STYLES["cell"]) for c in r] for r in rows]
            t = Table(data, hAlign="LEFT")
            t.setStyle(TableStyle([
                ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#4472C4")),
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#DCE6F1")),
                ("VALIGN", (0, 0), (-1, -1), "TOP")]))
            story += [Spacer(1, 4), t, Spacer(1, 8)]
        elif line.strip():
            para = line.strip()
            i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].startswith(("#", "|")):
                para += " " + lines[i].strip()
                i += 1
            i -= 1
            para = re.sub(r"\*([^*]+)\*", r"<i>\1</i>", para)
            story.append(Paragraph(para, STYLES["body"]))
        i += 1

    story.append(PageBreak())
    story.append(Paragraph("Figures", STYLES["h1"]))
    for fname, caption in FIGURES:
        path = os.path.join(ROOT, "figures", fname)
        if not os.path.exists(path):
            continue
        story.append(Image(path, width=5.6 * inch, height=5.6 * inch * 0.66))
        story.append(Paragraph(caption, STYLES["caption"]))

    story.insert(0, Paragraph(
        "Resistance-aware computational design of CRISPR gene-drive target sites in the "
        "<i>Anopheles gambiae</i> doublesex locus", STYLES["title"]))
    story.insert(1, Paragraph("Udita Phookan - 2026-09-24 - MEGA-27 program, item 15 (draft)", STYLES["author"]))
    story.insert(2, Spacer(1, 10))

    out = os.path.join(ROOT, "paper", "paper_draft.pdf")
    doc = SimpleDocTemplate(out, pagesize=letter, leftMargin=inch, rightMargin=inch,
                            topMargin=inch, bottomMargin=inch,
                            title="Resistance-aware computational design of CRISPR gene-drive target sites in dsx",
                            author="Udita Phookan")
    doc.build(story)
    return out


if __name__ == "__main__":
    print(build())
