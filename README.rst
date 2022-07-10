MolDrug
=======

|logo|

.. list-table::
    :widths: 12 35

    * - **Documentation**
      - |docs|
    * - **CI**
      - |tests|
    * - **Build**
      - |pypi-version|
    * - **Source Code**
      - |github|
    * - **Python Versions**
      - |pyversions|
    * - **Dependencies**
      - |rdkit| |crem| |openbabel|
    * - **License**
      - |license|
    * - **Downloads**
      - |downloads|


Description
-----------

**moldrug** is a python package for lead generation and optimization of small molecules. It use a Genetic Algorithm (GA) as searching engine in the chemical space and CReM library `crem <https://github.com/DrrDom/crem>`__ as chemical structure generator.

Documentation
-------------

The installation instructions, documentation and tutorials can be found online on `ReadTheDocs <https://moldrug.readthedocs.io/en/latest/>`_.

Issues
------

If you have found a bug, please open an issue on the `GitHub Issues <https://github.com/ale94mleon/moldrug/issues>`_.

Discussion
----------

If you have questions on how to use ProLIF, or if you want to give feedback or share ideas and new features, please head to the `GitHub Discussions <https://github.com/ale94mleon/moldrug/discussions>`_.

Citing **moldrug**
------------------

Please refer to the `citation page <https://moldrug.readthedocs.io/en/latest/source/citation.html>`__ on the documentation.

..  |logo|  image:: https://github.com/ale94mleon/moldrug/blob/main/row_data/logo.png?raw=true
    :target: https://github.com/ale94mleon/moldrug/
    :alt: logo
..  |docs|  image:: https://readthedocs.org/projects/moldrug/badge/?version=latest
    :target: https://moldrug.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation
..  |tests| image:: https://github.com/ale94mleon/moldrug/actions/workflows/conda.yml/badge.svg
    :target: https://github.com/ale94mleon/moldrug/actions/workflows/conda.yml/badge.svg
    :alt: CI/CD
..  |pypi-version|  image:: https://img.shields.io/pypi/v/moldrug.svg
    :target: https://pypi.python.org/pypi/moldrug/
    :alt: pypi-version
..  |github|    image:: https://badgen.net/badge/icon/github?icon=github&label
    :target: https://github.com/ale94mleon/moldrug
    :alt: GitHub-Repo
..  |pyversions|    image:: https://img.shields.io/pypi/pyversions/moldrug.svg
    :target: https://pypi.python.org/pypi/moldrug/
..  |rdkit| image:: https://img.shields.io/static/v1?label=Powered%20by&message=RDKit&color=3838ff&style=flat&logo=data:image/x-icon;base64,AAABAAEAEBAQAAAAAABoAwAAFgAAACgAAAAQAAAAIAAAAAEAGAAAAAAAAAMAABILAAASCwAAAAAAAAAAAADc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/FBT/FBT/FBT/FBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/PBT/PBT/PBT/PBT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/jIz/jIz/jIz/jIz/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/jIz/jIz/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/tLT/tLT/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/tLT/tLT/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/jIz/jIz/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/jIz/jIz/jIz/jIz/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/PBT/PBT/PBT/PBT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/FBT/FBT/FBT/FBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/////+B////AP///gB///wAP//4AB//+AAf//gAH//4AB//+AAf//gAH//8AD///gB///8A////gf////////
    :target: https://www.rdkit.org/docs/index.html
    :alt: rdkit
..  |openbabel| image:: https://img.shields.io/static/v1?label=Powered%20by&message=OpenBabel&color=6858ff&style=flat
    :target: https://openbabel.org/docs/dev/index.html
    :alt: openbabel
..  |crem| image:: https://img.shields.io/static/v1?label=Powered%20by&message=CReM&color=9438ff&style=flat
    :target: https://crem.readthedocs.io/en/latest/
    :alt: crem
..  |license| image:: https://badgen.net/pypi/license/moldrug/
    :target: https://pypi.python.org/pypi/moldrug/
    :alt: license
..  |downloads| image:: https://static.pepy.tech/personalized-badge/moldrug?period=month&units=international_system&left_color=grey&right_color=brightgreen&left_text=Downloads
    :target: https://pepy.tech/project/moldrug
    :alt: download
