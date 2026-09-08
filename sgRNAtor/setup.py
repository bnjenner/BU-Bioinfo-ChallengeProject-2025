import shutil
import warnings
from setuptools import setup, find_packages

# ── Check for required system binaries (bwa, samtools, bbmap) ──────────────
_REQUIRED_BINARIES = {
    "bwa":      "Install via conda  (bioconda::bwa)       or  https://github.com/lh3/bwa",
    "samtools": "Install via conda  (bioconda::samtools)  or  https://www.htslib.org",
    "bbduk.sh": "Install via conda  (bioconda::bbmap)     or  https://sourceforge.net/projects/bbmap/",
    "htstream": "Install via conda  (bioconda::htstream)  or  https://github.com/s4hts/HTStream"
}

_missing = {bin_: hint for bin_, hint in _REQUIRED_BINARIES.items()
            if shutil.which(bin_) is None}

if _missing:
    msg_lines = [
        "\n[sgRNAtor] The following required system tools were NOT found on PATH:"
    ]
    for bin_, hint in _missing.items():
        msg_lines.append(f"  • {bin_:<12}  →  {hint}")
    msg_lines.append(
        "\nThe Python package will be installed, but sgRNAtor will not run "
        "until these tools are available.\n"
        "Recommended: use the provided environment.yml with conda:\n"
        "  conda env create -f environment.yml && conda activate sgRNAtor\n"
    )
    warnings.warn("\n".join(msg_lines), stacklevel=2)

# ── Package metadata ───────────────────────────────────────────────────────
setup(
    name="sgRNAtor",
    version="0.0.1",
    description="Identification and quantification pipeline for sgRNA",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="B. N. Jenner",
    python_requires=">=3.10",

    # ── Pure-Python dependencies (installable by pip) ──────────────────────
    install_requires=[
        "biopython>=1.85",
        "editdistance>=0.8.1",
        "numpy>=2.0",
        "pysam>=0.22",          # compiles C extensions; needs htslib headers
    ],

    # ── Package discovery ──────────────────────────────────────────────────
    package_dir={"": "src"},
    packages=find_packages(where="src"),

    # ── CLI entry point ────────────────────────────────────────────────────
    entry_points={
        "console_scripts": [
            "sgRNAtor=sgRNAtor.main:main",
        ]
    },

    # ── Classifiers (optional but useful) ─────────────────────────────────
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)