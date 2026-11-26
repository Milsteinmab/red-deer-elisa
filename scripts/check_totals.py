"""Check the validation data against the values published in Hermida et al. (2020).

Standard library only. Run from the project root:

    python3 scripts/check_totals.py
"""

import csv
import math
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data" / "elisa_red_deer_validation.csv"

# Published values: Hermida et al. (2020), Small Ruminant Research 190, 106184, Table 2.
PUBLISHED = {
    "od_inhouse_elisa": {"cutoff": 0.26, "auc": 0.921, "se": 85.71, "sp": 88.89, "kappa": 0.62},
    "od_elisac": {"cutoff": 0.17, "auc": 0.603, "se": 57.14, "sp": 77.78, "kappa": 0.35},
}


def split(rows, column):
    pos = [float(r[column]) for r in rows if r["fecal_culture"] == "positive"]
    neg = [float(r[column]) for r in rows if r["fecal_culture"] == "negative"]
    return pos, neg


def table(pos, neg, cutoff, inclusive):
    above = (lambda x: x >= cutoff) if inclusive else (lambda x: x > cutoff)
    tp = sum(map(above, pos))
    fp = sum(map(above, neg))
    return tp, len(pos) - tp, fp, len(neg) - fp


def kappa(tp, fn, fp, tn):
    n = tp + fn + fp + tn
    po = (tp + tn) / n
    pe = ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n**2
    k = (po - pe) / (1 - pe)
    se = math.sqrt(po * (1 - po) / (n * (1 - pe) ** 2))
    return k, se


def auc(pos, neg):
    pairs = sum((p > q) + 0.5 * (p == q) for p in pos for q in neg)
    return pairs / (len(pos) * len(neg))


def main():
    with open(DATA, newline="") as f:
        rows = list(csv.DictReader(f))
    n_pos = sum(r["fecal_culture"] == "positive" for r in rows)
    n_neg = sum(r["fecal_culture"] == "negative" for r in rows)
    print(f"Samples: {len(rows)} ({n_pos} culture-positive, {n_neg} culture-negative)")
    assert (n_pos, n_neg) == (7, 9), "Expected 7 positive and 9 negative samples"

    for column, ref in PUBLISHED.items():
        pos, neg = split(rows, column)
        a = auc(pos, neg)
        print(f"\n{column}: AUC {a:.3f} (published {ref['auc']})")
        assert round(a, 3) == ref["auc"], "AUC does not match the published value"
        for inclusive in (False, True):
            tp, fn, fp, tn = table(pos, neg, ref["cutoff"], inclusive)
            k, se = kappa(tp, fn, fp, tn)
            sign = ">=" if inclusive else ">"
            print(
                f"  OD {sign} {ref['cutoff']}: TP {tp}, FN {fn}, FP {fp}, TN {tn} | "
                f"Se {100 * tp / (tp + fn):.2f}%, Sp {100 * tn / (tn + fp):.2f}%, "
                f"kappa {k:.3f} ± {se:.3f}"
            )
        print(
            f"  Published: Se {ref['se']}%, Sp {ref['sp']}%, kappa {ref['kappa']}"
        )


if __name__ == "__main__":
    main()
