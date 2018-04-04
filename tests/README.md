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
