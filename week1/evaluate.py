import subprocess
import time
import re
import os

DATASETS = ["data1", "data2", "data3", "data4"]

PYTHON_CMD = ["python", "week1/code/main.py"]
CODON_CMD = ["codon", "run", "-release", "week1/code/main.codon"]


def run_and_capture(cmd, dataset, language):
    """Run a command on a dataset, measure runtime, return (runtime_str, n50)."""
    start = time.time()
    try:
        result = subprocess.run(
            cmd + [f"week1/data/{dataset}"],
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running {language} on {dataset}")
        return "ERROR", 0
    except FileNotFoundError as e:
        print(f"❌ Command not found: {cmd[0]}")
        return "ERROR", 0

    end = time.time()
    elapsed = end - start

    # Format runtime as H:MM:SS
    h = int(elapsed // 3600)
    m = int((elapsed % 3600) // 60)
    s = int(elapsed % 60)
    runtime_str = f"{h}:{m:02d}:{s:02d}"

    # Try to find "N50 = <number>" in stdout
    n50 = 0
    match = re.search(r"N50\s*=\s*(\d+)", result.stdout)
    if match:
        n50 = int(match.group(1))

    return runtime_str, n50


def main():
    print("Dataset\tLanguage\tRuntime\tN50")
    print("-" * 80)

    for dataset in DATASETS:
        # Python version
        py_runtime, py_n50 = run_and_capture(PYTHON_CMD, dataset, "python")
        print(f"{dataset}\tpython\t{py_runtime}\t{py_n50}")

        # Codon version
        codon_runtime, codon_n50 = run_and_capture(CODON_CMD, dataset, "codon")
        print(f"{dataset}\tcodon\t{codon_runtime}\t{codon_n50}")


if __name__ == "__main__":
    main()
