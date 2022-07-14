# Package for Conda

Update version, URL in `meta.yaml`.

```
conda install conda-build
anaconda login
conda build -c conda-forge -c bioconda moldrug
```
Download the version and then:
sha256sum MolDrug-0.0.1.beta7.tar.gz
I should automate this thing

Valid `anaconda upload` command, with correct paths, is printed in the terminal at the end of building.

## For the documentation

```bash
make html
sphinx-build -b html source/. public
```

## Docker image

```bash
docker build -t 4moldrug . && docker run -it 4moldrug
docker images
docker stop $(docker ps -a -q)
docker rm $(docker ps -a -q)
docker rmi $(docker images -a -q)
docker logout
docker login
 ```
