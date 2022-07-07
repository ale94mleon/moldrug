FROM ubuntu:22.04

# https://serverfault.com/a/797318
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y
RUN apt-get upgrade -y

RUN apt-get install -y wget git

RUN useradd -m user

USER user
WORKDIR /home/user

RUN wget https://repo.anaconda.com/archive/Anaconda3-2021.04-Linux-x86_64.sh
RUN bash Anaconda3-2021.04-Linux-x86_64.sh -b -p ~/.anaconda

RUN echo "source ~/.anaconda/bin/activate" >> ~/.bashrc

SHELL ["/bin/bash", "-c"]

RUN source ~/.anaconda/bin/activate && conda install -y -c conda-forge rdkit
RUN source ~/.anaconda/bin/activate && conda install -y -c conda-forge openbabel
RUN source ~/.anaconda/bin/activate && conda install -y -c bioconda autodock-vina
