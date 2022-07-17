Installation
------------

Requirements:

    * Python 3.8+.
    * `RDKit <https://www.rdkit.org/docs/>`_ (2020.03+).
    * `Pandas <https://pandas.pydata.org/>`_.
    * `NumPy <https://numpy.org/>`_.
    * `tqdm <https://tqdm.github.io/>`_.
    * `CReM <https://github.com/DrrDom/crem>`_ (0.2.9+).
    * `Meeko <https://pypi.org/project/meeko/>`_.
    * `Vina <https://vina.scripps.edu/>`_.

.. note::

    If you have RDKit and Vina already installed you could try with ``pip install moldrug`` directlly.
    But if it is not the case or some version conflicts occurred, think about installed in a isoleated enviroment
    as it will be show in brief.

Via conda
~~~~~~~~~

It is recomendable to install through ``conda``:

.. code-block:: bash

    conda create -n moldrug
    conda activate moldrug
    conda install -c ale94mleon -c conda-forge -c bioconda moldrug

.. warning::

    Ussually pip has the lates stable version. But we are working to constantlly update the conda packege.
    Future plans are deployed inside conda-forge.

Via pip
~~~~~~~~~

Another possible way is direclly install from pip. But in this case you must have a correct installation
of RDKit and autodock-vina. One posibility is:

.. code-block:: bash

    conda create -n moldrug
    conda activate moldrug

Then install the dependencies libraries:

.. code-block:: bash

    conda install -y -c conda-forge rdkit">=2022.0"
    conda install -y -c bioconda autodock-vina

In the future we will consider to use the python modules `vina on pypi <https://pypi.org/project/vina/>`_. Finally:

.. code-block:: bash

    # To get the version on developing (not recomended)
    pip install git+https://github.com/ale94mleon/moldrug.git@main

or:

.. code-block:: bash

    # To get the last "stable" version (strongly recommended). This project is still in beta state.
    pip install moldrug

Work with a docker container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#. Use the `Docker configuration file on GitHub <https://github.com/ale94mleon/MolDrug/blob/main/Dockerfile>`__. 
#. Vist the `MolDrug <https://hub.docker.com/r/ale94mleon/4moldrug>`__ docker container.

Finaly ``pip install moldrug`` inside it.
