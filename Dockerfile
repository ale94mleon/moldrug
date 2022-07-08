FROM ubuntu:22.04

# https://serverfault.com/a/797318
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y
RUN apt-get upgrade -y

RUN apt-get install -y wget git sudo python3-pip

RUN echo "user ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers

RUN useradd -m user

USER user
WORKDIR /home/user

RUN wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
RUN bash Anaconda3-2022.05-Linux-x86_64.sh -b -p ~/.anaconda

SHELL ["/bin/bash", "-c"]
RUN source ~/.anaconda/bin/activate && conda init bash
RUN echo "export PATH=\$PATH:/home/user/.anaconda/bin" >> ~/.bashrc && source ~/.bashrc
RUN source ~/.anaconda/bin/activate && conda create -n druglead -y && echo "conda activate druglead" >> ~/.bashrc

RUN source ~/.anaconda/bin/activate && conda activate druglead && conda install -y -c conda-forge rdkit">=2022.0"
RUN source ~/.anaconda/bin/activate && conda activate druglead && conda install -y -c conda-forge openbabel">=3.1.0"
RUN source ~/.anaconda/bin/activate && conda activate druglead && conda install -y -c bioconda autodock-vina

# When druglead is release
# RUN source ~/.anaconda/bin/activate && conda activate druglead && pip install druglead