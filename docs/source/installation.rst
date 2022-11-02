Installation
============

Requirements:

    * `Python 3.8+ <https://docs.python.org/3/>`_.
    * `RDKit <https://www.rdkit.org/docs/>`_ (2020.03+).
    * `Pandas <https://pandas.pydata.org/>`_.
    * `NumPy <https://numpy.org/>`_.
    * `tqdm <https://tqdm.github.io/>`_.
    * `CReM <https://github.com/DrrDom/crem>`_ (0.2.9+).
    * `Meeko <https://pypi.org/project/meeko/>`_.
    * `Vina <https://vina.scripps.edu/>`_.

.. note::

    If you have RDKit and Vina already installed you could try with ``pip install moldrug`` directly.
    But if it is not the case or some version conflicts occurred, think about installed in a isolated environment
    as it will be show in brief.

Via pip (standard)
------------------

In this case you must have a correct installation
of RDKit and autodock-vina. If you already have it, skip the creation of the conda environment and the conda dependencies installation:

Create conda environment and install conda dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    conda create -n moldrug
    conda activate moldrug

Then install the dependencies libraries:

.. code-block:: bash

    conda install -y -c conda-forge rdkit">=2022.0" vina

..  In the future we will consider to use the python modules `vina on pypi <https://pypi.org/project/vina/>`_. Finally:

pip install
~~~~~~~~~~~

.. code-block:: bash

    # To get the last "stable" version (strongly recommended). This project is still in beta state.
    pip install moldrug

or:

.. code-block:: bash

    # To get the version on development (not recommended)
    pip install git+https://github.com/ale94mleon/MolDrug.git@main

Via conda
---------

We will create a new environment ``conda``:

.. code-block:: bash

    conda create -n moldrug
    conda activate moldrug
    conda install -c ale94mleon -c conda-forge moldrug

If some dependencies are missing, please installed through pip. Some of them could be:

.. code-block:: bash

    pip install meeko crem pyyaml scipy tqdm

.. note::

   In the future it will be deployed inside conda-forge.


Work with a docker container
----------------------------

#. Use the `Docker configuration file on GitHub <https://github.com/ale94mleon/MolDrug/blob/main/Dockerfile>`__.
#. Vist the `MolDrug <https://hub.docker.com/r/ale94mleon/4moldrug>`__ docker container.

Finally ``pip install moldrug`` inside it.
