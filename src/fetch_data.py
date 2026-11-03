"""Fetch public GDC data for Phase 1; use only standard-library modules.

Case-level expression is used only when exactly one RNA file exists for the case.
For cases with several files, extract the preselected PRIMARY TUMOR file instead.
This prevents normal tissue or repeated aliquots entering the patient-level row.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw"
CONFIG = json.loads((ROOT / "config/phase1.json").read_text())
API = CONFIG["api_base"]


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")


def download(endpoint, destination, payload=None):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size:
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    command = ["curl", "--silent", "--show-error", "--fail-with-body", "--retry", "5", "--retry-all-errors", "--retry-delay", "2",
               "--connect-timeout", "20", "--max-time", "180"]
    if payload is not None:
        request_path = RAW / "requests" / (destination.name + ".request.json")
        write_json(request_path, payload)
        command += ["-H", "Content-Type: application/json", "--data-binary", "@" + str(request_path)]
    command += [API + endpoint, "--output", str(temporary)]
    subprocess.run(command, check=True)
    temporary.replace(destination)


def load_hits(path):
    response = json.loads(path.read_text())
    warnings = response.get("warnings", {})
    if warnings:
        raise ValueError(f"GDC query warning: {warnings}")
    data = response["data"]
    if len(data["hits"]) != data["pagination"]["total"]:
        raise ValueError("GDC response is truncated; increase page size or paginate.")
    return data["hits"]


def sample_records(file):
    return [(case, sample) for case in file.get("cases", [])
            for sample in case.get("samples", [])]


def selection_key(file):
    samples = [s for _, s in sample_records(file) if s.get("sample_type") == CONFIG["sample_type"]]
    sample_id = min(s["submitter_id"] for s in samples)
    aliquots = [a.get("submitter_id", "") for s in samples
                for p in s.get("portions", []) for an in p.get("analytes", [])
                for a in an.get("aliquots", [])]
    return sample_id, min(aliquots, default=""), file["file_id"]


def select_files(files):
    by_case = defaultdict(list)
    for file in files:
        cases = {c["case_id"] for c in file.get("cases", [])}
        if len(cases) != 1:
            raise ValueError("A file has zero/multiple cases; inspect its mapping.")
        by_case[next(iter(cases))].append(file)
    selected = []
    for case_id, candidates in sorted(by_case.items()):
        primary = [f for f in candidates if any(s.get("sample_type") == CONFIG["sample_type"]
                                               for _, s in sample_records(f))]
        if not primary:
            continue
        chosen = min(primary, key=selection_key)
        selected.append({"case_id": case_id, "file_id": chosen["file_id"],
                         "file_name": chosen["file_name"], "md5sum": chosen["md5sum"],
                         "file_size": chosen["file_size"],
                         "sample_submitter_id": selection_key(chosen)[0],
                         "aliquot_submitter_id": selection_key(chosen)[1],
                         "all_rna_file_count": len(candidates),
                         "primary_rna_file_count": len(primary),
                         "expression_source": "selected_file" if len(candidates) > 1 else "single_file_case_api"})
    return selected


def direct_file(record):
    path = RAW / "rna" / (record["file_id"] + ".tsv")
    download("/data/" + record["file_id"], path)
    actual = hashlib.md5(path.read_bytes()).hexdigest()
    if actual != record["md5sum"]:
        raise ValueError(f"MD5 mismatch for {record['file_id']}; remove cached file and retry.")
    return record["file_id"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    download("/status", RAW / "gdc_status.json")
    cases_request = {
        "filters": {"op": "=", "content": {"field": "project.project_id", "value": CONFIG["project_id"]}},
        "size": 2000, "expand": "demographic,diagnoses",
        "fields": "case_id,submitter_id,project.project_id,demographic.demographic_id,demographic.sex_at_birth,demographic.race,demographic.ethnicity,demographic.vital_status,diagnoses.age_at_diagnosis,diagnoses.ajcc_pathologic_stage,diagnoses.diagnosis_id,diagnoses.classification_of_tumor,diagnoses.primary_diagnosis,diagnoses.tissue_or_organ_of_origin,diagnoses.submitter_id,diagnoses.days_to_diagnosis",
        "format": "JSON"}
    clauses = [("cases.project.project_id", CONFIG["project_id"]),
               ("data_type", "Gene Expression Quantification"),
               ("analysis.workflow_type", CONFIG["workflow"]), ("access", "open")]
    files_request = {
        "filters": {"op": "and", "content": [{"op": "=", "content": {"field": k, "value": v}} for k, v in clauses]},
        "size": 3000,
        "fields": "file_id,file_name,file_size,md5sum,state,data_type,data_category,experimental_strategy,access,analysis.workflow_type,cases.case_id,cases.submitter_id,cases.samples.sample_id,cases.samples.submitter_id,cases.samples.sample_type,cases.samples.tissue_type,cases.samples.tumor_descriptor,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id",
        "format": "JSON"}
    download("/cases", RAW / "clinical.json", cases_request)
    download("/files", RAW / "rna_files.json", files_request)
    cases = load_hits(RAW / "clinical.json")
    files = load_hits(RAW / "rna_files.json")
    selected = select_files(files)
    write_json(RAW / "selected_files.json", selected)
    print(f"Metadata: {len(cases)} cases, {len(files)} RNA files, {len(selected)} selected primary tumors.", flush=True)
    if args.metadata_only:
        return
    expression_request = {"case_ids": [r["case_id"] for r in selected],
                          "gene_ids": list(CONFIG["genes"].values()), "tsv_units": "uqfpkm", "format": "tsv"}
    download("/gene_expression/values", RAW / "case_expression.tsv", expression_request)
    # Deterministic validation of the API values for six cases with one RNA file.
    singleton_checks = [r for r in selected if r["all_rna_file_count"] == 1][:6]
    direct = [r for r in selected if r["all_rna_file_count"] > 1]
    needed = direct + singleton_checks
    print(f"Downloading {len(direct)} ambiguous-case selected files and {len(singleton_checks)} validation files.", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for n, result in enumerate(pool.map(direct_file, needed), 1):
            if n % 10 == 0 or n == len(needed):
                print(f"Validated RNA files: {n}/{len(needed)}", flush=True)
    write_json(RAW / "download_provenance.json", {
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "gdc_status": json.loads((RAW / "gdc_status.json").read_text()),
        "clinical_cases": len(cases), "rna_files": len(files),
        "selected_primary_tumor_cases": len(selected), "direct_files": len(direct),
        "singleton_validation_cases": [r["case_id"] for r in singleton_checks],
        "expression_unit": CONFIG["expression_unit"], "config": CONFIG})
    print("Public data ready. Run: python3 src/analyze_phase1.py", flush=True)


if __name__ == "__main__":
    main()
