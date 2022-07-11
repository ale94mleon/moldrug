# Package for Conda

Update version, URL in `meta.yaml`.

```
conda install conda-build
anaconda login
conda build -c conda-forge -c bioconda moldrug
```

Valid `anaconda upload` command, with correct paths, is printed in the terminal at the end of building.

## For the documentation

```bash
make html
sphinx-build -b html source/. public
```
