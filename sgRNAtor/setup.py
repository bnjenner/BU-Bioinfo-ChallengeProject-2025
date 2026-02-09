from setuptools import setup, find_packages

setup(
    name="sgRNAtor",
    version="0.0.1",
    description="Implements an sgRNA quantification pipeline",
    author="B. N. Jenner",
    python_requires=">=3.10",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    entry_points={
    "console_scripts": [
        "sgRNAtor=sgRNAtor.main:main",
    ]
}
)
