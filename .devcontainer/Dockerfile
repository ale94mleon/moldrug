FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       git wget curl ca-certificates sudo build-essential pandoc openssh-client \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN useradd -m -s /bin/bash user \
    && echo "user ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/user \
    && chmod 0440 /etc/sudoers.d/user

USER user
WORKDIR /home/user

# moldrug supports Python >=3.10,<3.12
ARG PYTHON_VERSION=3.10

# Miniforge (conda-forge by default, much lighter than Anaconda)
RUN wget -q "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh" -O miniforge.sh \
    && bash miniforge.sh -b -p /home/user/miniforge3 \
    && rm miniforge.sh

# The moldrug env first on PATH: works without `conda activate` (VS Code tasks, postCreateCommand)
ENV PATH="/home/user/miniforge3/envs/moldrug/bin:/home/user/miniforge3/bin:/home/user/.local/bin:${PATH}"

# Python is pinned explicitly; vina is the non-pip dependency of moldrug
RUN conda create -n moldrug -y -c conda-forge python=${PYTHON_VERSION} "rdkit>=2022.3.5" vina \
    && conda clean -afy

# Interactive shells: same env, with conda activated
SHELL ["/bin/bash", "-c"]
RUN conda init bash && echo "conda activate moldrug" >> ~/.bashrc

# Claude Code (native installer, installs into /home/user/.local/bin)
RUN curl -fsSL https://claude.ai/install.sh | bash

# Created here so the login volume mounted on it is owned by 'user'
RUN mkdir -p /home/user/.claude