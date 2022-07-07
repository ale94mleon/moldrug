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

RUN wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
RUN bash Anaconda3-2021.05-Linux-x86_64.sh -b -p ~/.anaconda

RUN echo "source ~/.anaconda/bin/activate" >> ~/.bashrc

SHELL ["/bin/bash", "-c"]

RUN source ~/.anaconda/bin/activate && conda create -n lead python=3.10

RUN source ~/.anaconda/bin/activate && conda activate lead && conda install -y -c conda-forge rdkit=2022.03.2
RUN source ~/.anaconda/bin/activate && conda activate lead && conda install -y -c conda-forge openbabel
RUN source ~/.anaconda/bin/activate && conda activate lead && conda install -y -c bioconda autodock-vina
