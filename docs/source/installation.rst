Installation
============

Requirements:

    * `Python 3.8+ <https://docs.python.org/3/>`_.
    * `RDKit <https://www.rdkit.org/docs/>`_ (2022.3.5+).
    * `Pandas <https://pandas.pydata.org/>`_.
    * `NumPy <https://numpy.org/>`_.
    * `tqdm <https://tqdm.github.io/>`_.
    * `CReM <https://github.com/DrrDom/crem>`_ (0.2.9+).
    * `Meeko <https://pypi.org/project/meeko/>`_.
    * `AutoDock-Vina <https://vina.scripps.edu/>`_.

.. note::
    `AutoDock-Vina <https://vina.scripps.edu/>`_ is the only non-pip dependency required for ``moldrug``. If you already have Vina installed, you can simply run the command ``pip install moldrug`` directly.
    However, in cases where Vina is not installed or if you encounter version conflicts, it's advisable to set up an isolated environment, either using `conda <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands>`_ or `python <https://docs.python.org/3/library/venv.html>`_.
    This ensures a clean and controlled environment for installing and running ``moldrug``.


Via pip (standard)
------------------
pip install
~~~~~~~~~~~

.. code-block:: bash

    # To get the last "stable" version (strongly recommended). This project is still in beta state.
    pip install moldrug

or:

.. code-block:: bash

    # To get the version on development (not recommended)
    pip install -U git+https://github.com/ale94mleon/MolDrug.git@main

There are multiple methods to obtain `AutoDock-Vina <https://vina.scripps.edu/>`_. You can use `conda <https://anaconda.org/conda-forge/vina>`_ or download the latest release from the `Vina Repository <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_.
We highly recommend using the ast release posted on method `Vina Repository <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_, which we will demonstrate here.


Getting last Vina
~~~~~~~~~~~~~~~~~

AutoDock-Vina is an ongoing project, and it is advisable to stay up-to-date by regularly checking for the latest `release <https://github.com/ccsb-scripps/AutoDock-Vina/releases/>`_.
As of the creation of this documentation, the most recent version is `v1.2.5 <https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.5>`_.
We will be using this version as a demonstration, but we strongly recommend using the most recent release for optimal performance and features. For detail information, please visit:
`AutoDock-Vina installation Instruction <https://autodock-vina.readthedocs.io/en/latest/installation.html>`_.

Unix- and Linux-based OS
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
    chmod a+x vina_1.2.5_linux_x86_64
    ./vina_1.2.5_linux_x86_64 -h

MacOS
^^^^^

.. code-block:: bash

    wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_mac_x86_64
    chmod a+x vina_1.2.5_mac_x86_64
    ./vina_1.2.5_mac_x86_64 -h

.. note::
    MacOs users might need to allow the execution of the application on ``Privacy & Security`` depending on the MacOS version.

Windows
^^^^^^^

Please, download form `release <https://github.com/ccsb-scripps/AutoDock-Vina/releases/>`_. Conda installation may not work.

Via conda
---------

MolDrug is also available through conda. However, the pip installation is the recommended one.

.. code-block:: bash

    conda create -n moldrug
    conda activate moldrug
    conda install -c ale94mleon -c conda-forge moldrug

.. note::
    MacOS users may face some problems trying to install because of the AutoDock-Vina dependency. If that is so, please follow the pip instructions.

If some dependencies are missing, please installed through pip. Some of them could be:

.. code-block:: bash

    pip install meeko crem pyyaml scipy tqdm


Work with a docker container
----------------------------

#. Use the `Docker configuration file on GitHub <https://github.com/ale94mleon/MolDrug/blob/main/Dockerfile>`__.
#. Vist the `MolDrug <https://hub.docker.com/r/ale94mleon/4moldrug>`__ docker container.

Finally ``pip install moldrug`` inside it.
