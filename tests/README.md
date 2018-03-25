# How to test

## C module

### Preparation

We need `googletest`.
We don't include googletest in this repo, so we need to download it from github.
The following does this,
```bash
./install_googletest.sh
```

### Run unittests

```bash
make test
```

### Prepare unittests
1. Make a test in `c` directory.
2. Add corresponding lines to Makefile.in

## Python module

### Preparation

We need `pytest` and `pytest-xdist`.
Install them via
```bash
conda install pytest pytest-xdist
```
for anaconda distribution or
```bash
pip install pytest pytest-xdist
```
for python native environment.

Run `run_pytest.sh`
```bash
./tests/run_pytest.sh
```
