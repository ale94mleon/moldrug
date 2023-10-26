MolDrug Dashboard
=================
|streamlit| |MolDrug| |rdkit| |mols2grid| |ProLIF| |py3Dmol| |stmol| |mdanalysis| |meeko| |numpy| |pandas| |seaborn|

`MolDrug-Dashboard <https://moldrug-dashboard.streamlit.app//>`__ will help you to get a overview 
of your MolDrug simulation.


Inputs
------
#. ``pbz2`` file. Exported when the command line is used. This file can be generated calling the function ``moldrug.utils.compressed_pickle``. The class ``moldrug.utils.Local`` and ``moldrug.utils.GA`` has already this method implemented.   
#. ``pdb`` file. The protein PDB if you would like to check the Protein-Ligand interaction network.

Home
----
The app presents a side bar (on the left) where the filters and representation options are located. The outputs are presented
in the center of the app.
|01|

Upload pbz2
-----------
The first step is to upload the pbz2 file. As soon as it is done some new options pops up in the side bar and the table of molecules is presented.
|02|

Customizing filters and properties
----------------------------------
All the properties used during the MolDrug run will be shown in ``Choose properties`` and can be selected.
The slide bar filter could be used to only show specific molecules in the table.
|03|

Interacting with the molecule table
-----------------------------------
Play around with the options of the table. You can:

* Display the properties of the molecule (depending on the ``Choose properties``).
* Get more properties making click on the picture.
* Sort based on other property (by default ``cost``).
* Highlight substructure with the SMART filter.

|04| 
|05|

Ligand-protein network interaction. Upload pdb
----------------------------------------------
To access this feature you must upload the pdb file. You could interact with the ProLIF image.
To change the molecule, simply introduce ``idx`` (the number at the top of each molecule the picture in the table)
of the desired molecule. 
|06|

Get 3D view
-----------
You can change the default representation to a 3D vie selecting ``py3Dmol`` in ``Representation``
|07|

Running info
------------
On this tab you can find valuable information about the convergency of MolDrug. It is customizable based on:
``Choose properties`` and ``Every how many generations`` bottoms.
|08|



.. Screenshots

..  |01| image:: _static/dashboard/01.png
    :alt: 01
..  |02| image:: _static/dashboard/02.png
    :alt: 02
..  |03| image:: _static/dashboard/03.png
    :alt: 03
..  |04| image:: _static/dashboard/04.png
    :alt: 04
..  |05| image:: _static/dashboard/05.png
    :alt: 05
..  |06| image:: _static/dashboard/06.png
    :alt: 06
..  |07| image:: _static/dashboard/07.png
    :alt: 07
..  |08| image:: _static/dashboard/08.png
    :alt: 08

.. Dependencies

..  |streamlit| image:: https://img.shields.io/static/v1?label=Powered%20by&message=Streamlit&color=FF3333&style=flat
    :target: https://streamlit.io/
    :alt: streamlit
..  |MolDrug| image:: https://img.shields.io/static/v1?label=Powered%20by&message=MolDrug&color=33FF5E&style=flat
    :target: https://moldrug.readthedocs.io/en/latest/
    :alt: MolDrug
..  |rdkit| image:: https://img.shields.io/static/v1?label=Powered%20by&message=RDKit&color=3838ff&style=flat&logo=data:image/x-icon;base64,AAABAAEAEBAQAAAAAABoAwAAFgAAACgAAAAQAAAAIAAAAAEAGAAAAAAAAAMAABILAAASCwAAAAAAAAAAAADc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/FBT/FBT/FBT/FBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/PBT/PBT/PBT/PBT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/jIz/jIz/jIz/jIz/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/jIz/jIz/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/tLT/tLT/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/tLT/tLT/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/jIz/jIz/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/jIz/jIz/jIz/jIz/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/PBT/PBT/PBT/PBT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/FBT/FBT/FBT/FBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/////+B////AP///gB///wAP//4AB//+AAf//gAH//4AB//+AAf//gAH//8AD///gB///8A////gf////////
    :target: https://www.rdkit.org/docs/index.html
    :alt: rdkit
..  |mols2grid| image:: https://img.shields.io/static/v1?label=Powered%20by&message=mols2grid&color=6858ff&style=flat
    :target: https://mols2grid.readthedocs.io/en/latest/
    :alt: mols2grid
..  |ProLIF| image:: https://img.shields.io/static/v1?label=Powered%20by&message=ProLIF&color=E933FF&style=flat
    :target:https://prolif.readthedocs.io/en/latest/
    :alt: ProLIF
..  |py3Dmol| image:: https://img.shields.io/static/v1?label=Powered%20by&message=py3Dmol&color=2BC3C5&style=flat
    :target:https://pypi.org/project/py3Dmol/
    :alt: py3Dmol
..  |stmol| image:: https://img.shields.io/static/v1?label=Powered%20by&message=stmol&color=D0E42B&style=flat
    :target:https://pypi.org/project/stmol/
    :alt: stmol
.. |mdanalysis| image:: https://img.shields.io/badge/Powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
    :alt: Powered by MDAnalysis
    :target: https://www.mdanalysis.org
..  |meeko| image:: https://img.shields.io/static/v1?label=Powered%20by&message=Meeko&color=56AD60&style=flat
    :target: https://github.com/forlilab/Meeko
    :alt: Meeko
..  |numpy| image:: https://img.shields.io/static/v1?label=Powered%20by&message=NumPy&color=49E1ED&style=flat
    :target: https://numpy.org/
    :alt: numpy
..  |pandas| image:: https://img.shields.io/static/v1?label=Powered%20by&message=Pandas&color=B2D7DA&style=flat
    :target: https://pandas.pydata.org/
    :alt: pandas
..  |seaborn|  image:: https://img.shields.io/static/v1?label=Powered%20by&message=seaborn&color=E523F5&style=flat
    :target: https://seaborn.pydata.org/
    :alt: seaborn