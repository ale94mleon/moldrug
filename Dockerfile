FROM conda/miniconda3

RUN conda update -n base -c defaults conda
RUN conda install -y -c conda-forge openbabel
RUN conda install -y -c bioconda autodock-vina
RUN conda install -y -c conda-forge rdkit cython
