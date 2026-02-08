from setuptools import setup, find_packages

setup(
    name="sgRNAtor",
    version="0.0.1",
    description="Implements an sgRNA quantification pipeline",
    author="B. N. Jenner",
    python_requires=">=3.10",
    packages=find_packages(include=["sgRNAtor", "sgRNAtor.*"]),
    entry_points={
        "console_scripts": [
            "sgRNAtor = sgRNAtor.sgRNAtor:main",
        ],
    },
)
