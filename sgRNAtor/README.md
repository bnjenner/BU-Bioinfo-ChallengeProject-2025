# sgRNAtor

## Installation
1. Install conda.

2. Create conda environment.
```
mkdir -p ./build
conda env create -f environment.yml --prefix $(pwd)/build
```

3. Activate your created conda environment

4. Install sgRNAtor by running this command in sgRNAtor directory containing setup.py.
```
pip install .
```

5. Test
```
sgRNAtor --help
```

