"""Patient-level, descriptive-only analysis of TCGA-BRCA Phase 1."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from fetch_data import CONFIG, ROOT, RAW, load_hits, select_files

os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "tmp/matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(ROOT / "tmp/cache"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = ROOT / "reports/phase1"
FIG = ROOT / "figures/phase1"
PROCESSED = ROOT / "data/processed"
GENES = list(CONFIG["genes"])
REVERSE_GENES = {v: k for k, v in CONFIG["genes"].items()}
STAGES = ["I", "II", "III", "IV"]


def write_csv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def numeric(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def summarize(values):
    a = np.asarray([v for v in values if v is not None and math.isfinite(v)], dtype=float)
    if not len(a):
        return {k: None for k in ["n", "mean", "median", "sd", "q1", "q3", "iqr", "min", "max", "outliers_iqr"]} | {"n": 0, "outliers_iqr": 0}
    q1, q3 = np.quantile(a, [.25, .75], method="linear")
    iqr = q3 - q1
    return {"n": len(a), "mean": float(a.mean()), "median": float(np.median(a)),
            "sd": float(a.std(ddof=1)) if len(a) > 1 else None,
            "q1": float(q1), "q3": float(q3), "iqr": float(iqr),
            "min": float(a.min()), "max": float(a.max()),
            "outliers_iqr": int(((a < q1 - 1.5 * iqr) | (a > q3 + 1.5 * iqr)).sum())}


def raw_expression(file_id):
    path = RAW / "rna" / (file_id + ".tsv")
    if not path.exists():
        raise FileNotFoundError(f"Missing selected file {file_id}; run fetch_data.py to completion.")
    result = {}
    with path.open() as stream:
        reader = csv.DictReader((line for line in stream if not line.startswith("#")), delimiter="\t")
        for row in reader:
            gene_id = row["gene_id"].split(".")[0]
            if gene_id in REVERSE_GENES:
                gene = REVERSE_GENES[gene_id]
                if row["gene_name"] != gene:
                    raise ValueError(f"Gene ID/symbol disagreement: {gene_id}, {row['gene_name']}, {gene}")
                if gene in result:
                    raise ValueError(f"Duplicate target gene in file {file_id}")
                result[gene] = numeric(row["fpkm_uq_unstranded"])
    if set(result) != set(GENES):
        raise ValueError(f"Selected genes missing from file {file_id}")
    return result


def case_expression():
    with (RAW / "case_expression.tsv").open() as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    case_ids = rows[0][1:]
    if len(case_ids) != len(set(case_ids)):
        raise ValueError("Duplicate case IDs in expression API result")
    values = {case_id: {} for case_id in case_ids}
    for row in rows[1:]:
        if len(row) != len(rows[0]):
            raise ValueError("Expression table is not rectangular")
        gene = REVERSE_GENES[row[0]]
        for case_id, value in zip(case_ids, row[1:]):
            values[case_id][gene] = numeric(value)
    return values


def main():
    for directory in [OUT, FIG, PROCESSED]:
        directory.mkdir(parents=True, exist_ok=True)
    cases = load_hits(RAW / "clinical.json")
    files = load_hits(RAW / "rna_files.json")
    selected = select_files(files)
    if len({c["case_id"] for c in cases}) != len(cases):
        raise ValueError("Duplicate cases in clinical response")
    if len({c["submitter_id"] for c in cases}) != len(cases):
        raise ValueError("Duplicate clinical patient barcodes")
    chosen = {r["case_id"]: r for r in selected}
    api_values = case_expression()
    if set(api_values) != set(chosen):
        raise ValueError("Expression API cases do not match selected molecular cases")

    validation = []
    for r in [r for r in selected if r["all_rna_file_count"] == 1][:6]:
        raw = raw_expression(r["file_id"])
        for gene in GENES:
            a, b = raw[gene], api_values[r["case_id"]].get(gene)
            ok = a is not None and b is not None and math.isclose(a, b, rel_tol=1e-7, abs_tol=1e-5)
            validation.append({"case_id": r["case_id"], "file_id": r["file_id"], "gene": gene,
                               "raw_fpkm_uq": a, "case_api_fpkm_uq": b, "passed": ok})
            if not ok:
                raise ValueError(f"Singleton expression validation failed: {r['case_id']} {gene}: {a} vs {b}")
    write_csv(OUT / "expression_validation.csv", validation)

    rows, diagnosis_audit, issues, sample_audit = [], [], [], []
    known_ids = {c["case_id"] for c in cases}
    if not set(chosen) <= known_ids:
        raise ValueError("Molecular case ID has no project case record")
    for c in sorted(cases, key=lambda c: c["submitter_id"]):
        case_id, patient = c["case_id"], c["submitter_id"]
        diagnoses = c.get("diagnoses", [])
        primary = [d for d in diagnoses if d.get("diagnosis_is_primary_disease") is True]
        diagnosis = primary[0] if len(primary) == 1 else {}
        has_clinical = bool(diagnoses or c.get("demographic"))
        diagnosis_audit.append({"case_id": case_id, "patient_id": patient,
                                "diagnosis_records": len(diagnoses), "primary_diagnosis_records": len(primary),
                                "selected_diagnosis_id": diagnosis.get("diagnosis_id", ""),
                                "selection_status": "unique_primary" if len(primary) == 1 else "missing_or_ambiguous_primary"})
        age_days = numeric(diagnosis.get("age_at_diagnosis"))
        age = age_days / 365.25 if age_days is not None else None
        if age is not None and not 0 <= age <= 120:
            issues.append({"patient_id": patient, "variable": "age_years", "value": age,
                           "decision": "invalid_age_set_missing"})
            age = None
        stage_raw = diagnosis.get("ajcc_pathologic_stage")
        stage_match = re.fullmatch(r"Stage (IV|III|II|I)(?:[ABC])?", str(stage_raw), flags=re.I)
        stage = stage_match.group(1).upper() if stage_match else None
        if stage is None and stage_raw not in [None, "", "Not Reported", "Unknown", "Not Applicable", "Stage X"]:
            issues.append({"patient_id": patient, "variable": "stage", "value": stage_raw,
                           "decision": "unrecognized_stage_set_missing"})
        r = chosen.get(case_id)
        values = (raw_expression(r["file_id"]) if r["all_rna_file_count"] > 1 else api_values[case_id]) if r else {}
        if r and not r["sample_submitter_id"].startswith(patient + "-"):
            raise ValueError(f"Patient/sample barcode mismatch: {patient} {r['sample_submitter_id']}")
        record = {"case_id": case_id, "patient_id": patient, "has_clinical": has_clinical,
                  "has_primary_rna": bool(r), "has_both": has_clinical and bool(r),
                  "age_years": age, "age_days_original": age_days,
                  "stage_original": stage_raw, "stage": stage,
                  "sex_at_birth": c.get("demographic", {}).get("sex_at_birth"),
                  "year_of_diagnosis": diagnosis.get("year_of_diagnosis"),
                  "age_is_obfuscated": c.get("demographic", {}).get("age_is_obfuscated"),
                  "diagnosis_id": diagnosis.get("diagnosis_id"),
                  "sample_id": r["sample_submitter_id"] if r else None,
                  "file_id": r["file_id"] if r else None,
                  "expression_source": r["expression_source"] if r else None}
        for gene in GENES:
            v = values.get(gene)
            if v is not None and (not math.isfinite(v) or v < 0):
                issues.append({"patient_id": patient, "variable": gene, "value": str(v),
                               "decision": "invalid_expression_set_missing"})
                v = None
            record[gene + "_fpkm_uq"] = v
            record[gene + "_log2"] = math.log2(v + 1) if v is not None else None
        record["complete_main_variables"] = bool(record["has_both"] and age is not None and stage is not None
                                                and all(record[g + "_log2"] is not None for g in GENES))
        rows.append(record)
    for f in files:
        for c in f.get("cases", []):
            for s in c.get("samples", []):
                selected_file = chosen.get(c["case_id"], {}).get("file_id") == f["file_id"]
                reason = "selected_primary_tumor" if selected_file else (
                    "non_primary_tissue" if s.get("sample_type") != "Primary Tumor" else "additional_primary_file")
                sample_audit.append({"case_id": c["case_id"], "patient_id": c["submitter_id"],
                                     "sample_id": s["sample_id"], "sample_barcode": s["submitter_id"],
                                     "sample_type": s.get("sample_type"), "file_id": f["file_id"],
                                     "file_name": f["file_name"], "source_md5": f["md5sum"],
                                     "source_bytes": f["file_size"], "selected": selected_file, "decision": reason})
    cohort = [r for r in rows if r["has_both"]]
    if not cohort:
        raise ValueError("Empty multimodal cohort")
    if len({r["case_id"] for r in cohort}) != len(cohort):
        raise ValueError("Patient-level cohort is not unique")
    write_csv(PROCESSED / "patients_phase1.csv", rows)
    write_csv(PROCESSED / "multimodal_cohort_phase1.csv", cohort)
    write_csv(OUT / "sample_file_audit.csv", sample_audit)
    write_csv(OUT / "diagnosis_audit.csv", diagnosis_audit)
    write_csv(OUT / "quality_issues.csv", issues, ["patient_id", "variable", "value", "decision"])
    write_csv(OUT / "selected_files.csv", selected)
    clinical = [r for r in rows if r["has_clinical"]]
    molecular = [r for r in rows if r["has_primary_rna"]]
    clinical_only = [r for r in rows if r["has_clinical"] and not r["has_primary_rna"]]
    molecular_only = [r for r in rows if r["has_primary_rna"] and not r["has_clinical"]]
    neither = [r for r in rows if not r["has_primary_rna"] and not r["has_clinical"]]
    complete = [r for r in cohort if r["complete_main_variables"]]
    flow = [
        {"step": "project_case_records", "n": len(rows)},
        {"step": "cases_with_clinical_records", "n": len(clinical)},
        {"step": "clinical_diagnosis_records_before_selection", "n": sum(len(c.get("diagnoses", [])) for c in cases)},
        {"step": "cases_with_unique_primary_diagnosis", "n": sum(r["primary_diagnosis_records"] == 1 for r in diagnosis_audit)},
        {"step": "open_star_counts_files_all_tissues", "n": len(files)},
        {"step": "primary_tumor_files_before_patient_selection", "n": sum(s["sample_type"] == "Primary Tumor" for s in sample_audit)},
        {"step": "non_primary_files_excluded", "n": sum(s["sample_type"] != "Primary Tumor" for s in sample_audit)},
        {"step": "additional_primary_files_excluded", "n": sum(s["decision"] == "additional_primary_file" for s in sample_audit)},
        {"step": "unique_primary_rna_patients", "n": len(molecular)},
        {"step": "clinical_only", "n": len(clinical_only)},
        {"step": "molecular_only", "n": len(molecular_only)},
        {"step": "neither_modality", "n": len(neither)},
        {"step": "main_multimodal_patient_cohort", "n": len(cohort)},
        {"step": "complete_all_main_variables", "n": len(complete)}]
    write_csv(OUT / "data_flow.csv", flow)
    stats = []
    quantitative = ["age_years"] + [g + suffix for g in GENES for suffix in ["_fpkm_uq", "_log2"]]
    for population, group, variables in [("clinical_before_merge", clinical, ["age_years"]),
                                          ("molecular_before_merge", molecular, quantitative[1:]),
                                          ("multimodal_cohort", cohort, quantitative)]:
        for variable in variables:
            summary = summarize([r[variable] for r in group])
            stats.append({"population": population, "variable": variable, "population_n": len(group),
                          "missing_n": len(group) - summary["n"],
                          "missing_percent": 100 * (len(group) - summary["n"]) / len(group), **summary})
    write_csv(OUT / "descriptive_statistics.csv", stats)
    categorical = []
    for population, group in [("clinical_before_merge", clinical), ("multimodal_cohort", cohort)]:
        for variable in ["stage_original", "stage", "sex_at_birth"]:
            counts = Counter(r[variable] or "Missing" for r in group)
            for category, n in sorted(counts.items()):
                categorical.append({"population": population, "variable": variable, "category": category,
                                    "n": n, "denominator": len(group), "percent": 100 * n / len(group)})
    write_csv(OUT / "categorical_frequencies.csv", categorical)
    grouped = []
    for stage in STAGES:
        group = [r for r in cohort if r["stage"] == stage]
        grouped.append({"stage": stage, "group_n": len(group),
                        **summarize([r["MKI67_log2"] for r in group])})
    write_csv(OUT / "mki67_by_stage.csv", grouped)
    availability = []
    for label, group in [("both", cohort), ("clinical_only", clinical_only),
                         ("molecular_only", molecular_only), ("neither", neither)]:
        availability.append({"group": label, "n": len(group), **{"age_" + k: v for k, v in summarize([r["age_years"] for r in group]).items()}})
    write_csv(OUT / "modality_availability.csv", availability)
    selection_comparison = []
    for label, group in [("both", cohort), ("clinical_only", clinical_only)]:
        for variable in ["stage", "sex_at_birth"]:
            counts = Counter(r[variable] or "Missing" for r in group)
            for category in sorted(set(r[variable] or "Missing" for r in clinical)):
                n = counts[category]
                selection_comparison.append({"group": label, "variable": variable, "category": category,
                                             "n": n, "denominator": len(group),
                                             "percent": 100 * n / len(group) if group else None})
    write_csv(OUT / "clinical_selection_comparison.csv", selection_comparison)

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": "white", "axes.titleweight": "bold", "savefig.dpi": 180})
    colors = ["#256b83", "#d78b43", "#52865d", "#9a627e"]
    manifest = []

    def save(name, fig, title, caption):
        fig.tight_layout()
        fig.savefig(FIG / (name + ".png"), dpi=180, bbox_inches="tight")
        fig.savefig(FIG / (name + ".svg"), bbox_inches="tight")
        plt.close(fig)
        manifest.append({"file": f"figures/phase1/{name}.png", "title": title, "caption": caption})

    a = [r["age_years"] for r in cohort if r["age_years"] is not None]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.hist(a, bins=20, color=colors[0], edgecolor="white")
    ax.axvline(np.median(a), color=colors[1], lw=2, label=f"Median = {np.median(a):.1f} years")
    ax.set(xlabel="Age at primary diagnosis (years)", ylabel="Patients (count)", title=f"Age distribution | available n = {len(a)}")
    ax.legend(frameon=False)
    save("01_age_histogram", fig, "توزیع سن هنگام تشخیص", "سن از تشخیص اصلی استخراج و از روز به سال تبدیل شده است. خط عمودی میانه را نشان می‌دهد؛ داده‌های گمشده در هیستوگرام وارد نمی‌شوند.")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    stages = STAGES + ["Missing"]
    counts = [sum((r["stage"] or "Missing") == s for r in cohort) for s in stages]
    bars = ax.bar(stages, counts, color=colors[0])
    ax.bar_label(bars, padding=3)
    ax.set_ylim(0, max(counts) * 1.14)
    ax.set(xlabel="Pathological stage (ordered groups)", ylabel="Patients (count)", title=f"Disease stage | cohort n = {len(cohort)}")
    save("02_stage_counts", fig, "فراوانی مرحلهٔ بیماری", "زیرمرحله‌های معتبر در چهار مرحلهٔ اصلی ادغام شده‌اند. ستون Missing نیز نمایش داده شده تا مخرج درصدها و نبود داده روشن بماند.")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3))
    for ax, variable, title in zip(axes, ["MKI67_fpkm_uq", "MKI67_log2"], ["Original normalized scale", "Log-transformed scale"]):
        v = [r[variable] for r in cohort if r[variable] is not None]
        ax.hist(v, bins=30, color=colors[0], edgecolor="white")
        ax.set(title=title, xlabel="MKI67 (FPKM-UQ)" if variable.endswith("uq") else "MKI67: log2(FPKM-UQ + 1)", ylabel="Patients (count)")
    fig.suptitle("MKI67 distribution", fontweight="bold", y=1.03)
    save("03_mki67_distribution", fig, "توزیع بیان MKI67 در دو مقیاس", "نمایش هم‌زمان مقیاس اصلی و تبدیل لگاریتمی، اثر چولگی و مقادیر بزرگ را آشکار می‌کند. تبدیل لگاریتمی آزمون نرمال‌بودن یا تضمین نرمال‌بودن نیست.")
    fig, ax = plt.subplots(figsize=(8, 4.7))
    data = [[r["MKI67_log2"] for r in cohort if r["stage"] == s and r["MKI67_log2"] is not None] for s in STAGES]
    valid_stages = [(s, v) for s, v in zip(STAGES, data) if v]
    ax.boxplot([v for s, v in valid_stages], tick_labels=[f"{s}\nn={len(v)}" for s, v in valid_stages], patch_artist=True,
               boxprops={"facecolor": "#d6e8ee"}, medianprops={"color": "#bd6b29", "linewidth": 2},
               flierprops={"marker": ".", "markersize": 4})
    ax.set(xlabel="Pathological stage", ylabel="MKI67: log2(FPKM-UQ + 1)", title="MKI67 by stage | descriptive comparison")
    save("04_mki67_by_stage", fig, "مقایسهٔ توصیفی MKI67 میان مراحل بیماری", "خط وسط میانه، جعبه فاصلهٔ چارک اول تا سوم و نقاط بیرونی مقادیر خارج از حصار ۱٫۵ برابر IQR را نشان می‌دهند. این نقاط حذف نشده‌اند؛ اندازهٔ هر گروه زیر محور آمده است.")
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    for ax, gene, color in zip(axes.flat, GENES, colors):
        v = [r[gene + "_log2"] for r in cohort if r[gene + "_log2"] is not None]
        ax.hist(v, bins=25, color=color, edgecolor="white")
        ax.set(title=f"{gene} | n={len(v)}", xlabel="log2(FPKM-UQ + 1)", ylabel="Patients (count)")
    save("05_selected_genes", fig, "توزیع چهار ژن منتخب", "همهٔ ژن‌ها پیش از محاسبهٔ نتایج انتخاب شده‌اند. هر پنل توزیع یک ژن را در همان مقیاس مشخص نشان می‌دهد و هیچ آزمون انتخاب ویژگی اجرا نشده است.")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    group_names = ["Both modalities", "Clinical only", "Molecular only", "Neither"]
    group_values = [len(cohort), len(clinical_only), len(molecular_only), len(neither)]
    bars = ax.bar(group_names, group_values, color=[colors[0], colors[1], colors[2], "#9ca3af"])
    ax.bar_label(bars, padding=3)
    ax.set_ylim(0, max(group_values) * 1.13)
    ax.set(ylabel="Patients (count)", xlabel="Data availability", title=f"Modality availability | project cases n={len(rows)}")
    save("06_modality_availability", fig, "الگوی دسترسی به انواع داده", "دادهٔ بالینی به معنی وجود رکورد تشخیص یا جمعیت‌شناختی است؛ وجود رکورد به معنی کامل‌بودن همهٔ متغیرها نیست. دادهٔ مولکولی به نمونهٔ تومور اولیهٔ واجد شرایط اشاره دارد.")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    visible = ["age_years", "stage"] + [g + "_fpkm_uq" for g in GENES]
    missing = [100 * sum(r[v] is None for r in cohort) / len(cohort) for v in visible]
    bars = ax.barh(["Age", "Stage"] + GENES, missing, color=colors[0])
    ax.bar_label(bars, labels=[f"{v:.2f}%" for v in missing], padding=4)
    ax.set_xlim(0, max(10, max(missing) * 1.4))
    ax.set(xlabel="Missing observations (% of multimodal cohort)", title="Missingness in main variables")
    ax.invert_yaxis()
    save("07_missingness", fig, "درصد گمشدگی متغیرهای اصلی", "مخرج همهٔ درصدهای این شکل، تعداد بیماران دارای هر دو نوع داده است. گمشدگی متغیر با نبود کل یک نوع داده یکسان نیست.")
    dump(OUT / "figure_manifest.json", manifest)
    years = [r["year_of_diagnosis"] for r in rows if isinstance(r["year_of_diagnosis"], int)]
    status = json.loads((RAW / "gdc_status.json").read_text())
    fingerprints = []
    for path in sorted(RAW.rglob("*")):
        if path.is_file() and not path.name.endswith(".part"):
            fingerprints.append({"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size,
                                 "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    write_csv(OUT / "source_fingerprints.csv", fingerprints)
    provenance = json.loads((RAW / "download_provenance.json").read_text())
    dump(OUT / "provenance.json", provenance)
    results = {"generated_at_utc": datetime.now(timezone.utc).isoformat(), "gdc_status": status,
               "flow": flow, "sample_types": dict(Counter(r["sample_type"] for r in sample_audit)),
               "n_main": len(cohort), "n_complete": len(complete),
               "quantitative": stats, "categorical": categorical, "mki67_by_stage": grouped,
               "availability": availability, "selection_comparison": selection_comparison,
               "quality_issue_count": len(issues), "validation_checks": len(validation),
               "age_obfuscated_n": sum(r["age_is_obfuscated"] is True for r in cohort),
               "diagnosis_year_min": min(years) if years else None, "diagnosis_year_max": max(years) if years else None,
               "diagnosis_year_available_n": len(years),
               "primary_diagnosis_missing_or_ambiguous": sum(r["primary_diagnosis_records"] != 1 for r in diagnosis_audit),
               "all_rna_files_per_case": dict(Counter(r["all_rna_file_count"] for r in selected)),
               "final_patient_ids_unique": len({r["case_id"] for r in cohort}) == len(cohort),
               "phase_scope": "descriptive_only_phase_1", "figures": manifest}
    dump(OUT / "results.json", results)
    print(json.dumps({"project_cases": len(rows), "multimodal_patients": len(cohort),
                      "complete_all_main_variables": len(complete), "figures": len(manifest),
                      "validation_checks_passed": len(validation), "quality_issues": len(issues)}, indent=2))


if __name__ == "__main__":
    main()
