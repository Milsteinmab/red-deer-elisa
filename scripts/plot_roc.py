"""Regenerate figures/elisa_roc.svg with Python's standard library.

From the project root: python3 scripts/plot_roc.py
"""

import csv
from collections import Counter
from pathlib import Path

from check_totals import DATA, PUBLISHED, split, table


OUTPUT = Path(__file__).resolve().parent.parent / "figures" / "elisa_roc.svg"


def roc_points(positive, negative):
    """Move all sera with the same OD together when lowering the threshold."""
    cases = Counter(positive)
    controls = Counter(negative)
    tp = fp = 0
    points = [(0.0, 0.0)]
    for score in sorted(cases.keys() | controls.keys(), reverse=True):
        tp += cases[score]
        fp += controls[score]
        points.append((fp / len(negative), tp / len(positive)))
    return points


def curve_area(points):
    return sum(
        (x2 - x1) * (y1 + y2) / 2
        for (x1, y1), (x2, y2) in zip(points, points[1:])
    )


def plot_xy(fpr, tpr):
    # Square plotting area; SVG's vertical coordinate increases downward.
    return 84 + 400 * fpr, 500 - 400 * tpr


def polyline(points, style):
    coords = " ".join(f"{x:.3f},{y:.3f}" for x, y in (plot_xy(*p) for p in points))
    return f'<polyline class="{style}" fill="none" points="{coords}"/>'


def render(curves, areas, counts, cutoff):
    tp, fn, fp, tn = counts
    sensitivity = tp / (tp + fn)
    specificity = tn / (fp + tn)
    cx, cy = plot_xy(fp / (fp + tn), sensitivity)
    ticks = []
    for tick in (0, 0.25, 0.5, 0.75, 1):
        x, y = plot_xy(tick, tick)
        ticks.extend([
            f'<path class="grid" d="M{x} 100 V500 M84 {y} H484"/>',
            f'<text x="{x}" y="523" text-anchor="middle">{tick:g}</text>',
            f'<text x="71" y="{y + 5}" text-anchor="end">{tick:g}</text>',
        ])

    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="800" height="610" viewBox="0 0 800 610" role="img" aria-labelledby="title desc">
  <title id="title">ELISA performance in red deer</title>
  <desc id="desc">Empirical ROC curves against fecal culture. In-house ELISA AUC {areas["od_inhouse_elisa"]:.3f}, bovine ELISA AUC {areas["od_elisac"]:.3f}. An open circle marks the in-house threshold OD greater than {cutoff}, with sensitivity {sensitivity:.1%} and specificity {specificity:.1%}.</desc>
  <style>
    text {{ font-family: Arial, Helvetica, sans-serif; fill: #172c38; font-size: 15px }}
    .heading {{ font-size: 28px; font-weight: bold }}
    .label {{ font-weight: bold }}
    .grid {{ fill: none; stroke: #dae1e5; stroke-width: 1 }}
    .inhouse {{ stroke: #005c83; stroke-width: 3.5 }}
    .bovine {{ stroke: #a04417; stroke-width: 3; stroke-dasharray: 9 5 }}
    .cutoff {{ fill: white; stroke: #005c83; stroke-width: 3 }}
  </style>
  <rect width="800" height="610" fill="white"/>
  <text class="heading" x="32" y="42">ELISA performance in red deer</text>
  <text x="32" y="70">Fecal culture reference · {tp + fn} positive and {tn + fp} negative sera</text>
  {"".join(ticks)}
  <path d="M84 100 V500 H484" fill="none" stroke="#516572"/>
  <path d="M84 500 L484 100" stroke="#85929b" stroke-dasharray="3 5"/>
  {polyline(curves["od_inhouse_elisa"], "inhouse")}
  {polyline(curves["od_elisac"], "bovine")}
  <circle class="cutoff" cx="{cx:.3f}" cy="{cy:.3f}" r="6.5"/>
  <text x="284" y="553" text-anchor="middle">False-positive rate (1 − specificity)</text>
  <text transform="translate(26 300) rotate(-90)" text-anchor="middle">Sensitivity</text>
  <path class="inhouse" d="M520 151 h38"/>
  <text class="label" x="570" y="156">In-house ELISA</text>
  <text x="570" y="179">AUC {areas["od_inhouse_elisa"]:.3f}</text>
  <path class="bovine" d="M520 219 h38"/>
  <text class="label" x="570" y="224">Bovine ELISA (ELISAc)</text>
  <text x="570" y="247">AUC {areas["od_elisac"]:.3f}</text>
  <circle class="cutoff" cx="539" cy="291" r="6.5"/>
  <text class="label" x="570" y="296">In-house OD &gt; {cutoff}</text>
  <text x="570" y="320">Sensitivity {tp}/{tp + fn} ({sensitivity:.1%})</text>
  <text x="570" y="344">Specificity {tn}/{fp + tn} ({specificity:.1%})</text>
  <text x="32" y="589">Hermida et al. (2020) · Small Ruminant Research 190, 106184</text>
</svg>
'''


def main():
    with DATA.open(newline="", encoding="utf-8") as source:
        rows = list(csv.DictReader(source))

    curves, areas = {}, {}
    for assay, reference in PUBLISHED.items():
        positive, negative = split(rows, assay)
        if (len(rows), len(positive), len(negative)) != (16, 7, 9):
            raise ValueError("Expected 7 culture-positive and 9 culture-negative sera")
        points = roc_points(positive, negative)
        area = curve_area(points)
        if round(area, 3) != reference["auc"]:
            raise ValueError(f"{assay} ROC area differs from the published AUC")
        curves[assay] = points
        areas[assay] = area

    positive, negative = split(rows, "od_inhouse_elisa")
    cutoff = PUBLISHED["od_inhouse_elisa"]["cutoff"]
    counts = table(positive, negative, cutoff, inclusive=False)
    if counts != (6, 1, 1, 8):
        raise ValueError("In-house counts differ from the verified strict-cutoff table")
    svg = render(curves, areas, counts, cutoff)
    OUTPUT.parent.mkdir(exist_ok=True)
    OUTPUT.write_text(svg, encoding="utf-8")
    print(f"Wrote {OUTPUT.relative_to(OUTPUT.parent.parent)}")


if __name__ == "__main__":
    main()
