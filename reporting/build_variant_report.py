#!/usr/bin/env python3
"""
Turns this pipeline's report_variants.tsv / report_summary.tsv into a single
polished HTML page (print it to PDF from the browser, Cmd/Ctrl+P) suitable
for presenting a run's results.

Deliberately does NOT try to reproduce sarek's own QC report: sarek already
writes a full MultiQC report at <sarek_outdir>/multiqc/multiqc_report.html
covering alignment rate, duplicates, coverage, etc. Pass --multiqc to link
it from this page instead of re-deriving those numbers here.

Usage:
  python3 build_variant_report.py \
      --report-tsv results/wholeexome/report_variants/<sample>/report/report_variants.tsv \
      --summary-tsv results/wholeexome/report_variants/<sample>/report/report_summary.tsv \
      --sample-id VIDO-027-01 \
      --genome GRCh38 \
      --exome-note "generic UCSC RefSeq all-exons target set (real capture kit not yet confirmed)" \
      --multiqc results/sarek/multiqc/multiqc_report.html \
      --out variant_report.html

Only needs the Python standard library, so it runs anywhere, no venv needed.
"""
import argparse
import csv
import html
import sys
from collections import Counter
from datetime import date
from pathlib import Path

ACCENT = "#2e4a52"
ACCENT_SOFT = "#4d6b73"

CURATED_COLUMNS = [
    ("POSITION", "Position"),
    ("VARIANT_GENE_NAME", "Gene"),
    ("VARIANT_ANNOTATION", "Consequence"),
    ("VARIANT_HGVS_C", "HGVS.c"),
    ("VARIANT_HGVS_P", "HGVS.p"),
    ("SAMPLE_GT", "Genotype"),
    ("SAMPLE_DEPTH", "Depth"),
    ("SAMPLE_ALT_AF", "Alt AF"),
    ("CLINVAR_SIG", "ClinVar"),
    ("GNOMAD_AF", "gnomAD AF"),
    ("CADD_PHRED", "CADD"),
    ("REVEL", "REVEL"),
]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def esc(x):
    return html.escape("" if x is None else str(x))


def bar_chart(counts, title, unit_label):
    """counts: list of (label, value), already sorted. Single-hue horizontal bars."""
    if not counts:
        return f'<p class="chart-empty">No data for "{title}".</p>'
    max_val = max(v for _, v in counts) or 1
    row_h = 26
    gap = 10
    label_w = 190
    bar_area = 360
    height = len(counts) * (row_h + gap) - gap + 30
    width = label_w + bar_area + 60
    rows = []
    y = 20
    for label, val in counts:
        bar_w = max(2, (val / max_val) * bar_area)
        rows.append(f'''
        <text x="{label_w - 8}" y="{y + row_h / 2 + 4}" text-anchor="end"
              class="chart-label">{esc(label)}</text>
        <rect x="{label_w}" y="{y}" width="{bar_w:.1f}" height="{row_h - 6}"
              rx="3" fill="{ACCENT}"/>
        <text x="{label_w + bar_w + 8}" y="{y + row_h / 2 + 4}" class="chart-value">{val}</text>
        ''')
        y += row_h + gap
    svg = f'''
    <svg viewBox="0 0 {width} {height}" width="100%" height="{height}"
         role="img" aria-label="{esc(title)}: {unit_label}">
      {''.join(rows)}
    </svg>
    '''
    return svg


def build(args):
    report_rows = read_tsv(args.report_tsv)
    summary_rows = read_tsv(args.summary_tsv) if args.summary_tsv else []

    sample_id = args.sample_id or (report_rows[0].get("SAMPLE_ID") if report_rows else "sample")

    n_variants = len(report_rows)

    ann_counts = Counter(r.get("VARIANT_ANNOTATION", "unknown") or "unknown" for r in report_rows)
    ann_sorted = sorted(ann_counts.items(), key=lambda kv: -kv[1])[:10]

    depth_bins = Counter()
    for r in report_rows:
        try:
            d = float(r.get("SAMPLE_DEPTH") or 0)
        except ValueError:
            d = 0
        if d < 10:
            b = "< 10x"
        elif d < 20:
            b = "10-19x"
        elif d < 30:
            b = "20-29x"
        elif d < 50:
            b = "30-49x"
        elif d < 100:
            b = "50-99x"
        else:
            b = "100x+"
        depth_bins[b] += 1
    bin_order = ["< 10x", "10-19x", "20-29x", "30-49x", "50-99x", "100x+"]
    depth_sorted = [(b, depth_bins[b]) for b in bin_order if depth_bins[b] > 0]

    table_rows = []
    for r in report_rows:
        cells = "".join(f"<td>{esc(r.get(col, ''))}</td>" for col, _ in CURATED_COLUMNS)
        table_rows.append(f"<tr>{cells}</tr>")

    summary_html = ""
    if summary_rows:
        headers = list(summary_rows[0].keys())
        head = "".join(f"<th>{esc(h)}</th>" for h in headers)
        body = "".join(
            "<tr>" + "".join(f"<td>{esc(r.get(h, ''))}</td>" for h in headers) + "</tr>"
            for r in summary_rows
        )
        summary_html = f'''
        <div class="table-wrap">
          <table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>
        </div>
        '''

    multiqc_note = ""
    if args.multiqc:
        multiqc_note = f'''
        <p>Alignment and coverage QC for this sample is in sarek's own
        MultiQC report, at <span class="mono">{esc(args.multiqc)}</span>.
        That report is not reproduced here since it already covers this
        thoroughly; open it alongside this page.</p>
        '''
    else:
        multiqc_note = '''
        <p>Alignment and coverage QC for this sample is in sarek's own
        MultiQC report (<span class="mono">&lt;sarek_outdir&gt;/multiqc/multiqc_report.html</span>),
        not reproduced here.</p>
        '''

    doc = f'''<title>{esc(sample_id)} Variant Report</title>
<style>
  :root {{
    --bg: #fdfcfa; --surface: #f7f5f0; --text: #1c1c1a; --text-muted: #55534d;
    --rule: #d9d5c9; --rule-strong: #b8b3a3; --accent: {ACCENT}; --accent-soft: {ACCENT_SOFT};
    --table-border: #a8a397;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      --bg: #17181a; --surface: #1f2123; --text: #e9e7e0; --text-muted: #b0ada2;
      --rule: #3a3c3d; --rule-strong: #52544f; --accent: #8fb7bf; --accent-soft: #7a9ea6;
      --table-border: #4d4f4a;
    }}
  }}
  :root[data-theme="dark"] {{
    --bg: #17181a; --surface: #1f2123; --text: #e9e7e0; --text-muted: #b0ada2;
    --rule: #3a3c3d; --rule-strong: #52544f; --accent: #8fb7bf; --accent-soft: #7a9ea6;
    --table-border: #4d4f4a;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    background: var(--bg); color: var(--text);
    font-family: Georgia, "Times New Roman", Times, serif;
    font-size: 16px; line-height: 1.6;
    padding-inline: 24px; padding-block: 48px 72px;
  }}
  .sheet {{ max-width: 800px; margin: 0 auto; }}
  .label {{ font-size: 12px; letter-spacing: 0.11em; text-transform: uppercase; color: var(--accent-soft); margin: 0 0 8px; }}
  h1 {{ font-size: 27px; margin: 0 0 6px; text-wrap: balance; }}
  .dateline {{ font-size: 14px; color: var(--text-muted); font-style: italic; margin: 0 0 4px; }}
  hr.rule {{ border: none; border-top: 1px solid var(--rule-strong); margin: 28px 0; }}
  h2 {{ font-size: 19px; margin: 0 0 12px; padding-bottom: 6px; border-bottom: 1px solid var(--rule); }}
  section {{ margin-bottom: 32px; }}
  p {{ margin: 0 0 13px; }}
  p:last-child {{ margin-bottom: 0; }}
  .mono {{ font-family: "Courier New", Courier, monospace; font-size: 0.9em; background: var(--surface); padding: 0.1em 0.35em; border-radius: 2px; }}
  .table-wrap {{ overflow-x: auto; margin: 0 0 8px; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 13.5px; }}
  th, td {{ border: 1px solid var(--table-border); padding: 6px 8px; text-align: left; white-space: nowrap; }}
  th {{ background: var(--surface); font-weight: 700; }}
  .chart-label {{ font-size: 12px; fill: var(--text); font-family: Georgia, serif; }}
  .chart-value {{ font-size: 12px; fill: var(--text-muted); font-family: Georgia, serif; }}
  .chart-empty {{ color: var(--text-muted); font-style: italic; }}
  .caveat {{ padding: 12px 16px; border-left: 2px solid var(--rule-strong); color: var(--text-muted); font-size: 14.5px; }}
  .figure {{ margin: 0 0 26px; }}
  .figure-title {{ font-weight: 700; font-size: 14.5px; margin: 0 0 8px; }}
</style>
<div class="sheet">
  <p class="label">Variant Report</p>
  <h1>{esc(sample_id)}: whole exome variant report</h1>
  <p class="dateline">Generated {date.today().isoformat()} &middot; nextflow-wholeexome &middot; genome build {esc(args.genome)}</p>

  <hr class="rule">

  <section>
    <h2>Run notes</h2>
    <p>This report covers {n_variants} variant{'s' if n_variants != 1 else ''} passing this
    pipeline's default filters (minimum depth 10, minimum alternate allele depth 3, minimum
    alternate allele frequency 5%), annotated with ClinVar, COSMIC, dbSNP, gnomAD, CADD, REVEL,
    and HPO gene-disease associations, following the same reporting logic our ngstk exome module
    has always used.</p>
    <p class="caveat">{esc(args.exome_note)}</p>
    {multiqc_note}
  </section>

  <hr class="rule">

  <section>
    <h2>Summary</h2>
    {summary_html or '<p class="chart-empty">No summary file provided.</p>'}
  </section>

  <hr class="rule">

  <section>
    <h2>Figures</h2>
    <div class="figure">
      <p class="figure-title">Variants by predicted consequence</p>
      {bar_chart(ann_sorted, "Variants by predicted consequence", "variant count")}
    </div>
    <div class="figure">
      <p class="figure-title">Variants by read depth at the called site</p>
      {bar_chart(depth_sorted, "Variants by read depth", "variant count")}
    </div>
  </section>

  <hr class="rule">

  <section>
    <h2>Variant table</h2>
    <div class="table-wrap">
      <table>
        <thead><tr>{"".join(f"<th>{esc(label)}</th>" for _, label in CURATED_COLUMNS)}</tr></thead>
        <tbody>{"".join(table_rows) or "<tr><td>No variants.</td></tr>"}</tbody>
      </table>
    </div>
    <p style="font-size: 13px; color: var(--text-muted); margin-top: 8px;">
      {len(CURATED_COLUMNS)} of {len(report_rows[0]) if report_rows else 0} columns shown here;
      the full report_variants.tsv has the rest (gnomAD subpopulation frequencies, SpliceAI,
      PrimateAI, dbscSNV, segmental duplication and repeat-region flags, HPO disease terms).
    </p>
  </section>
</div>
'''
    Path(args.out).write_text(doc)
    print(f"Wrote {args.out} ({n_variants} variants)", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--report-tsv", required=True)
    p.add_argument("--summary-tsv")
    p.add_argument("--sample-id")
    p.add_argument("--genome", default="GRCh38")
    p.add_argument("--exome-note",
                    default="Annotated against the generic UCSC RefSeq all-exons target set, "
                            "not yet a capture-kit-specific manifest.")
    p.add_argument("--multiqc")
    p.add_argument("--out", required=True)
    build(p.parse_args())


if __name__ == "__main__":
    main()
