"""Compute probability tables and figures from the archived TCGA-BRCA cohort."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports/phase2"
FIGURES = ROOT / "figures/phase2"
CLINICAL_GROUPS = ("I–II", "III–IV")
MARKER_GROUPS = ("At or below median", "Above median")
COLORS = ("#307296", "#be5828")


def numeric(value):
    return float(value) if value else np.nan


def table(name, header, rows):
    with (REPORTS / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def save(fig, name):
    for fmt in ("png", "svg"):
        fig.savefig(FIGURES / f"{name}.{fmt}", dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def hypergeometric(k, n, hits, size):
    if not 0 <= k <= hits or not 0 <= n - k <= size - hits:
        return 0.0
    return math.comb(hits, k) * math.comb(size - hits, n - k) / math.comb(size, n)


def main():
    config = json.loads((ROOT / "config/phase2.json").read_text())
    source = ROOT / config["input"]
    with source.open(encoding="utf-8-sig", newline="") as stream:
        cohort = list(csv.DictReader(stream))
    for key in ("patient_id", "case_id"):
        identifiers = [row[key] for row in cohort]
        if not all(identifiers) or len(set(identifiers)) != len(cohort):
            raise ValueError(f"Expected exactly one nonempty {key} per patient")
    if any(row["has_both"] != "True" or not row["sample_id"] for row in cohort):
        raise ValueError("Cohort must have linked clinical data and primary tumor RNA")

    marker = np.array([numeric(row[config["molecular_variable"]]) for row in cohort])
    age = np.array([numeric(row[config["clinical_variable"]]) for row in cohort])
    stage = np.array([row["stage"] for row in cohort])
    if set(stage) - {"", *config["early_stage"], *config["advanced_stage"]}:
        raise ValueError("Unexpected pathological stage")
    has_marker = np.isfinite(marker)
    has_age = np.isfinite(age)
    if np.any(marker[has_marker] < 0) or np.any((age[has_age] < 0) | (age[has_age] > 120)):
        raise ValueError("Invalid marker or age on the Phase 1 scales")
    x = marker[has_marker]
    if x.size < 2 or x.std(ddof=0) == 0:
        raise ValueError("Marker distribution must have at least two distinct observations")
    median = float(np.median(x))
    high_marker = marker > median
    valid_stage = stage != ""
    event_mask = has_marker & valid_stage
    age_mask = has_marker & has_age
    advanced = np.isin(stage, config["advanced_stage"])
    counts = np.array([
        [np.count_nonzero(event_mask & (advanced == clinical) & (high_marker == molecular))
         for molecular in (False, True)]
        for clinical in (False, True)
    ])
    n_event = int(np.sum(event_mask))
    if counts.sum() != n_event or min(*counts.sum(axis=0), *counts.sum(axis=1)) == 0 or age_mask.sum() < 2:
        raise ValueError("Empty event margin or insufficient paired observations")
    joint = counts / n_event
    p_advanced = float(joint[1, :].sum())
    p_high = float(joint[:, 1].sum())
    p_both = float(joint[1, 1])
    p_stage_given_high = float(counts[1, 1] / counts[:, 1].sum())
    p_stage_given_low = float(counts[1, 0] / counts[:, 0].sum())
    p_high_given_stage = float(counts[1, 1] / counts[1, :].sum())
    p_high_given_early = float(counts[0, 1] / counts[0, :].sum())
    bayes = p_high_given_stage * p_advanced / p_high
    if not math.isclose(bayes, p_stage_given_high):
        raise ValueError("Bayes calculation differs from direct counting")

    n_draws = config["binomial_draws"]
    if n_draws < 1 or n_draws > n_event:
        raise ValueError("Invalid count-model sample size")
    support = list(range(n_draws + 1))
    binom = [math.comb(n_draws, k) * p_high**k * (1 - p_high)**(n_draws - k) for k in support]
    exact = [hypergeometric(k, n_draws, int(counts[:, 1].sum()), n_event) for k in support]
    if not all(math.isclose(sum(pmf), 1) for pmf in (binom, exact)):
        raise ValueError("A probability mass function did not sum to one")

    mean = float(x.mean())
    sd = float(x.std(ddof=0))
    q1, q3 = np.quantile(x, [.25, .75])
    iqr = q3 - q1
    low_fence, high_fence = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    standardized = (x - mean) / sd
    normal_q = np.array([NormalDist().inv_cdf((i - .5) / x.size) for i in range(1, x.size + 1)])
    fit_q = mean + sd * normal_q
    empirical_q = np.sort(x)
    paired_age, paired_x = age[age_mask], marker[age_mask]
    paired_age_mean, paired_x_mean = paired_age.mean(), paired_x.mean()
    covariance = float(np.mean((paired_age - paired_age_mean) * (paired_x - paired_x_mean)))
    covariance_sample = covariance * paired_age.size / (paired_age.size - 1)
    results = {
        "source": config["input"], "release": config["input_release"],
        "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "cohort_n": len(cohort), "event_n": n_event, "marker_n": int(x.size), "pair_n": int(paired_x.size),
        "missing_stage": int((~valid_stage).sum()), "missing_marker": int((~has_marker).sum()),
        "missing_age": int((~has_age).sum()), "median_ties": int(np.count_nonzero(x == median)),
        "counts": counts.tolist(), "joint": joint.tolist(),
        "p_advanced": p_advanced, "p_high": p_high, "p_both": p_both,
        "p_stage_given_high": p_stage_given_high, "p_stage_given_low": p_stage_given_low,
        "p_high_given_stage": p_high_given_stage, "p_high_given_early": p_high_given_early,
        "independent_joint": p_advanced * p_high, "joint_gap": p_both - p_advanced * p_high,
        "bayes": bayes, "stage_given_high_minus_low": p_stage_given_high - p_stage_given_low,
        "n_draws": n_draws, "binomial_mean": n_draws * p_high,
        "binomial_var": n_draws * p_high * (1 - p_high),
        "hypergeometric_var": n_draws * p_high * (1 - p_high) * (n_event - n_draws) / (n_event - 1),
        "pmf_total_variation": .5 * sum(abs(a - b) for a, b in zip(binom, exact)),
        "median": median, "mean": mean, "variance_empirical": sd**2,
        "variance_sample": float(x.var(ddof=1)), "sd": sd, "q1": float(q1), "q3": float(q3),
        "min": float(x.min()), "max": float(x.max()),
        "skewness": float(np.mean(standardized**3)),
        "lower_fence": float(low_fence), "upper_fence": float(high_fence),
        "lower_outliers": int(np.sum(x < low_fence)), "upper_outliers": int(np.sum(x > high_fence)),
        "normal_below_zero": NormalDist(mean, sd).cdf(0),
        "age_mean_paired": float(paired_age_mean), "marker_mean_paired": float(paired_x_mean),
        "covariance_empirical": covariance, "covariance_sample": covariance_sample,
    }
    REPORTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    (REPORTS / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    table("event_table.csv", ("stage_group", "marker_group", "count", "joint_probability"),
          [(clinical, molecular, counts[i, j], joint[i, j])
           for i, clinical in enumerate(CLINICAL_GROUPS) for j, molecular in enumerate(MARKER_GROUPS)])
    table("count_pmf.csv", ("k", "binomial_with_replacement", "hypergeometric_without_replacement"),
          zip(support, binom, exact))

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(6.2, 3.4), layout="constrained")
    ax.bar(MARKER_GROUPS, [p_stage_given_low, p_stage_given_high], color=COLORS, width=.55)
    ax.axhline(p_advanced, color="0.35", linestyle="--", label="Overall within analysis cohort")
    for i, value in enumerate((p_stage_given_low, p_stage_given_high)):
        ax.text(i, value + .006, f"{value:.1%}", ha="center")
    ax.set(ylabel="Fraction at pathological stage III–IV", xlabel="MKI67 expression group",
           ylim=(0, min(1, max(p_stage_given_low, p_stage_given_high) * 1.3 + .02)))
    ax.legend(loc="upper left", frameon=False)
    save(fig, "01_conditional_probability")

    fig, ax = plt.subplots(figsize=(6.2, 3.4), layout="constrained")
    ax.bar(support, binom, color=COLORS[0], alpha=.85, label="With replacement (binomial)")
    ax.plot(support, exact, "o", color=COLORS[1], markersize=4, label="Without replacement (exact)")
    ax.set(xlabel=f"Patients above MKI67 median in {n_draws} draws", ylabel="Probability")
    ax.legend(frameon=False)
    save(fig, "02_count_model")

    fig, (left, right) = plt.subplots(1, 2, figsize=(8.4, 3.4), layout="constrained")
    left.hist(x, bins="fd", density=True, color=COLORS[0], alpha=.83)
    support_x = np.linspace(float(x.min()), float(x.max()), 250)
    density = np.exp(-.5 * ((support_x - mean) / sd)**2) / (sd * np.sqrt(2 * np.pi))
    left.plot(support_x, density, color=COLORS[1], label="Fitted normal")
    left.set(xlabel="MKI67 log₂(FPKM-UQ + 1)", ylabel="Density", title="Histogram")
    left.legend(frameon=False)
    right.scatter(fit_q, empirical_q, s=7, color=COLORS[0])
    lo, hi = min(fit_q[0], empirical_q[0]), max(fit_q[-1], empirical_q[-1])
    right.plot([lo, hi], [lo, hi], color=COLORS[1], linestyle="--")
    right.set(xlabel="Fitted normal quantiles", ylabel="Observed MKI67 quantiles", title="Normal Q–Q")
    save(fig, "03_mki67_distribution")

    fig, ax = plt.subplots(figsize=(6.2, 3.7), layout="constrained")
    ax.scatter(paired_age, paired_x, alpha=.35, s=11, color=COLORS[0], edgecolors="none")
    ax.axvline(paired_age_mean, color="0.5", linestyle="--", linewidth=.8)
    ax.axhline(paired_x_mean, color="0.5", linestyle="--", linewidth=.8)
    ax.set(xlabel="Age at diagnosis (years)", ylabel="MKI67 log₂(FPKM-UQ + 1)")
    save(fig, "04_age_mki67")
    print(f"Analyzed {len(cohort)} patients. Stage–marker n={n_event}; age–marker n={paired_x.size}.")


if __name__ == "__main__":
    main()
