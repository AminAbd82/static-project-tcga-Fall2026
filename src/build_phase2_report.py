"""Render the Persian Phase 2 report using XeLaTeX and the analysis results."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]


def main():
    results = json.loads((ROOT / "reports/phase2/results.json").read_text())
    xelatex = shutil.which("xelatex")
    if not xelatex:
        raise SystemExit("Install XeLaTeX and the xepersian package to build the report.")
    template = (ROOT / "reports/phase2/report_fa.tex").read_text()
    values = {key: f"{value:.4f}" if isinstance(value, float) else str(value)
              for key, value in results.items() if isinstance(value, (int, float))}
    for key, prefix in (("counts", "Count"), ("joint", "Joint")):
        for i, row in enumerate(results[key]):
            for j, value in enumerate(row):
                values[f"{prefix}{i}{j}"] = str(value) if key == "counts" else f"{100 * value:.2f}"
    for key in ("p_advanced", "p_high", "p_both", "p_stage_given_high",
                "p_stage_given_low", "p_high_given_stage", "p_high_given_early",
                "independent_joint", "joint_gap", "bayes",
                "stage_given_high_minus_low", "normal_below_zero"):
        values[f"pct_{key}"] = f"{100 * results[key]:.2f}"
    values.update(
        early_total=sum(results["counts"][0]), advanced_total=sum(results["counts"][1]),
        low_total=results["counts"][0][0] + results["counts"][1][0],
        high_total=results["counts"][0][1] + results["counts"][1][1],
        release=results["release"], input_sha256=results["input_sha256"],
        normal_direction=("تقریب نرمال برای مرکز توزیع مناسب است، اما دنباله‌ها را فقط تقریبی نشان می‌دهد."
                          if abs(results["skewness"]) < .5 else
                          "به علت چولگی، تقریب نرمال باید با احتیاط به کار رود."),
        covariance_direction=("منفی" if results["covariance_empirical"] < 0 else
                              "مثبت" if results["covariance_empirical"] > 0 else "صفر"),
    )
    for key, value in values.items():
        template = template.replace(f"@@{key}@@", str(value))
    if "@@" in template:
        raise ValueError("Unresolved result in the report")
    with TemporaryDirectory(prefix="tcga-phase2-") as tmp:
        tex = Path(tmp) / "phase2.tex"
        tex.write_text(template)
        env = dict(os.environ, openin_any="p", openout_any="p")
        # Use the repository root for relative font and figure paths.
        for _ in range(2):
            run = subprocess.run(
                [xelatex, "-interaction=nonstopmode", "-halt-on-error",
                 f"-output-directory={tmp}", str(tex)],
                cwd=ROOT, env=env, capture_output=True, text=True,
            )
            if run.returncode:
                raise SystemExit(run.stdout + run.stderr)
        shutil.copyfile(Path(tmp) / "phase2.pdf", ROOT / "reports/phase2.pdf")


if __name__ == "__main__":
    main()
