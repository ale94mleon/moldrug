# moldrug Dashboard

[![streamlit](https://img.shields.io/static/v1?label=Powered%20by&message=Streamlit&color=FF3333&style=flat)](https://streamlit.io/)
[![moldrug](https://img.shields.io/static/v1?label=Powered%20by&message=moldrug&color=33FF5E&style=flat)](https://moldrug.readthedocs.io/en/latest/)
[![rdkit](https://img.shields.io/static/v1?label=Powered%20by&message=RDKit&color=3838ff&style=flat&logo=data:image/x-icon;base64,AAABAAEAEBAQAAAAAABoAwAAFgAAACgAAAAQAAAAIAAAAAEAGAAAAAAAAAMAABILAAASCwAAAAAAAAAAAADc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/FBT/FBT/FBT/FBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/PBT/PBT/PBT/PBT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/jIz/jIz/jIz/jIz/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/jIz/jIz/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/tLT/tLT/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/tLT/tLT/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/jIz/jIz/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/jIz/jIz/jIz/jIz/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/PBT/PBT/PBT/PBT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/FBT/FBT/FBT/FBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/////+B////AP///gB///wAP//4AB//+AAf//gAH//4AB//+AAf//gAH//8AD///gB///8A////gf////////)](https://www.rdkit.org/docs/index.html)
[![mols2grid](https://img.shields.io/static/v1?label=Powered%20by&message=mols2grid&color=6858ff&style=flat)](https://mols2grid.readthedocs.io/en/latest/)
[![ProLIF](https://img.shields.io/static/v1?label=Powered%20by&message=ProLIF&color=E933FF&style=flat)](https://prolif.readthedocs.io/en/latest/)
[![py3Dmol](https://img.shields.io/static/v1?label=Powered%20by&message=py3Dmol&color=2BC3C5&style=flat)](https://pypi.org/project/py3Dmol/)
[![stmol](https://img.shields.io/static/v1?label=Powered%20by&message=stmol&color=D0E42B&style=flat)](https://pypi.org/project/stmol/)
[![mdanalysis](https://img.shields.io/badge/Powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)
[![meeko](https://img.shields.io/static/v1?label=Powered%20by&message=Meeko&color=56AD60&style=flat)](https://github.com/forlilab/Meeko)
[![numpy](https://img.shields.io/static/v1?label=Powered%20by&message=NumPy&color=49E1ED&style=flat)](https://numpy.org/)
[![pandas](https://img.shields.io/static/v1?label=Powered%20by&message=Pandas&color=B2D7DA&style=flat)](https://pandas.pydata.org/)
[![seaborn](https://img.shields.io/static/v1?label=Powered%20by&message=seaborn&color=E523F5&style=flat)](https://seaborn.pydata.org/)

## Getting ready

[moldrug-Dashboard](https://moldrug-dashboard.streamlit.app/) will help you to get an overview of your moldrug simulation. It can be run online, but it is preferable to run it locally to not run with memory issues. To do so, you have to install its dependencies (you could activate your moldrug environment or create a new one).

```bash
pip install -r https://raw.githubusercontent.com/ale94mleon/moldrug/main/streamlit/requirements.txt
```

And run the application with

```bash
streamlit run https://raw.githubusercontent.com/ale94mleon/moldrug/main/streamlit/moldrug-dashboard.py
```

Of course, you can always download [moldrug-dashboard.py](https://github.com/ale94mleon/moldrug/blob/main/streamlit/moldrug-dashboard.py) from the GitHub repository to your personal computer.

## Moldrug-Dashboard at a glance

![01](_static/dashboard/dashboard.svg)

## Inputs

- `pbz2` file. Exported when the command line is used. This file can be generated by calling the function {py:func}`moldrug.utils.compressed_pickle`. The classes {py:class}`moldrug.utils.Local` and {py:class}`moldrug.utils.GA` have already this method implemented (`pickle`).
- `pdb` file. The protein PDB if you would like to check the Protein-Ligand interaction network. In this case the `pdbqt` attribute of the {py:class}`moldrug.utils.Individual` shoud have been updated to the docking pose. This is done automatically by the classes {py:class}`moldrug.utils.GA` and {py:class}`moldrug.utils.Local`

## Upload pbz2

The first step is to upload the pbz2 file. As soon as it is done some new options pops up in the sidebar (**A**) and the table of molecules is presented.

## Customizing filters and properties

All the properties used during the moldrug run will be shown in ``Choose properties`` and can be selected. The slide bar filter could be used to only show specific molecules in the table.

## Interacting with the table of molecules

Play around with the options of the table. You can:

- Display the properties of the molecule (depending on the ``Choose properties``).
- Get more properties making click on the picture.
- Sort based on other property (by default ``cost``).
- Highlight substructure with the SMART filter.
- A summary table of properties, sortable and downloadable, is also accessible (**B**).

## Ligand-protein network interaction

To access this feature you must upload the pdb file (**C**). You can interact with the ProLIF image. To change the molecule, simply introduce `idx` (the number at the top of each molecule the picture in the table) of the desired molecule.

A comprehensive table of interaction types for all solutions is provided, enabling sorting, filtering, and download functionalities (**E**)

## Set 3D view

You can change the default 2D representation to a 3D view in `Representation` (**D**). This representation is also interactive.

## Running info

This tab (**F**) offers a view of the population's evolution across generations. This information is pivotal to assess convergence of the moldrug simulation. It is customizable based on: ``Choose properties`` and ``Every how many generations`` bottoms.
