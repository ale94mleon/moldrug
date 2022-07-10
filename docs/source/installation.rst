Installation
------------

Requirements:

* Python 3.8+
* `RDKit <https://www.rdkit.org/docs/>`_ (2020.03+)
* `Pandas <https://pandas.pydata.org/>`_ (1.0+)
* `NumPy <https://numpy.org/>`_
* `sklearn <https://scikit-learn.org/stable/>`_
* `tqdm <https://tqdm.github.io/>`_
* `CReM <https://github.com/DrrDom/crem>`_ (0.2.9+)
* `OpenBabel <https://openbabel.org/docs/dev/Installation/install.html>`_ (3.1.0+)

We are currently working in a ``conda`` and a ``docker`` container.

For now to get a fully functional MolDrug is recommendable the following::

    # create a separate virtual environment
    conda create -n moldrug
    # activate it
    conda activate moldrug

Then simply run the following command to install the dependencies libraries::

    conda install -y -c conda-forge rdkit">=2022.0"
    conda install -y -c conda-forge openbabel">=3.1.0"
    conda install -y -c bioconda autodock-vina

If you already have them; just skip the last steep. Finally::

    # To get the version on developing
    pip install git+https://github.com/ale94mleon/moldrug.git@main

or

    # To get the last "stable" version. This project is still in beta state.
    pip install moldrug
