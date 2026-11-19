# An image with moldrug and its only non-pip dependency, AutoDock-Vina, already installed.
#
# moldrug is a plain Python package, so `pip install moldrug` is usually enough. This image
# exists for the awkward part: AutoDock-Vina is not on PyPI, and installing it by hand is
# unpleasant on some systems, MacOS with Gatekeeper in particular.
#
# Build:
#   docker build -t moldrug .
#   docker build -t moldrug --build-arg MOLDRUG_VERSION=4.0.0 .   # pin a release
#
# Run:
#   docker run --rm -it -v "$PWD":/home/user/work moldrug
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
# moldrug supports Python >=3.10,<3.14
ARG PYTHON_VERSION=3.11
# Empty installs the latest release; set it to pin one
ARG MOLDRUG_VERSION=

# The versions are deliberately not pinned. Debian keeps only the current version of a
# package, so a pin breaks the build as soon as one is superseded, and for ca-certificates
# and curl the patched version is the one the image should carry.
# hadolint ignore=DL3008
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates curl bzip2 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN useradd -m -s /bin/bash user
USER user
WORKDIR /home/user

# Miniforge for the architecture the image is being built on, so this works on x86_64 and
# on arm64 alike. The previous version of this file hardcoded an x86_64 Anaconda installer
# and could not build on arm64 at all.
RUN curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh" -o miniforge.sh \
    && bash miniforge.sh -b -p /home/user/miniforge3 \
    && rm miniforge.sh

# The environment first on PATH, so it is used without `conda activate`
ENV PATH="/home/user/miniforge3/envs/moldrug/bin:/home/user/miniforge3/bin:${PATH}"

# vina is the point of this image; rdkit comes from conda-forge to avoid building it
RUN conda create -n moldrug -y -c conda-forge \
        python=${PYTHON_VERSION} rdkit vina \
    && conda clean -afy

RUN pip install --no-cache-dir "moldrug${MOLDRUG_VERSION:+==${MOLDRUG_VERSION}}"

# Fail the build rather than ship an image where either half is broken
RUN vina --version && python -c "import moldrug; print('moldrug', moldrug.__version__)"

# Interactive shells get the same environment, with conda activated
SHELL ["/bin/bash", "-c"]
RUN conda init bash && echo "conda activate moldrug" >> ~/.bashrc

WORKDIR /home/user/work
CMD ["bash"]
