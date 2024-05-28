import os
import tempfile
import time
import urllib.parse
from io import StringIO

import matplotlib.pyplot as plt
import MDAnalysis as mda
import mols2grid
import numpy as np
import pandas as pd
import prolif as plf
import pubchempy as pcp
import requests
import seaborn as sns
from bs4 import BeautifulSoup
from meeko import PDBQTMolecule, RDKitMolCreate
from pandas.api.types import (is_categorical_dtype, is_datetime64_any_dtype,
                              is_numeric_dtype, is_object_dtype)
from rdkit import Chem, DataStructs
from stmol import showmol

import streamlit as st
import streamlit.components.v1 as components
from moldrug import utils

# TODO
# add SyGma for metabolic prediction
# add prediction of synthetic routes.
# st.set_page_config('wide')
st.sidebar.empty()
st.title('Dashboard')
st.image('https://raw.githubusercontent.com/ale94mleon/MolDrug/main/docs/source/_static/MolDrug-logo-full.svg', width=200)

with st.expander('**About the App**'):
    st.markdown("👈 Open the side bar to introduce the data.\n\n"
                "This app is to get an overview of a **MolDrug**'s result at glance. "
                "Check this [flash tutorial](https://moldrug.readthedocs.io/en/latest/source/moldrug_dahsboard.html) "
                "in case you get stock on how to use the app; or [MolDrug's docs](https://moldrug.rtfd.io/) "
                "and [MolDrug's GitHub](https://github.com/ale94mleon/moldrug/) for more information.\n\n"
                "This project received foundings from [Marie Skłodowska-Curie Actions](https://cordis.europa.eu/project/id/860592). "
                "It was developed in the [Computational Biophysics Group](https://biophys.uni-saarland.de/) "
                "of [Saarland University](https://www.uni-saarland.de/en/home.html) in collaboration "
                "with the pharmaceutical company [Boehringer Ingelheim](https://www.boehringer-ingelheim.com/de/).")

tab1, tab2, tab3 = st.tabs(["Molecules", "Running info", "Novelty Checker"])


@st.cache_data
def convert_df(df):
    return df.to_csv().encode('utf-8')


def filter_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Take it from: https://blog.streamlit.io/auto-generate-a-dataframe-filtering-ui-in-streamlit-with-filter_dataframe/
    Adds a UI on top of a dataframe to let viewers filter columns

    Args:
        df (pd.DataFrame): Original dataframe

    Returns:
        pd.DataFrame: Filtered dataframe
    """
    modify = st.checkbox("Add filters")

    if not modify:
        return df

    df = df.copy()

    # Try to convert datetimes into a standard format (datetime, no timezone)
    for col in df.columns:
        if is_object_dtype(df[col]):
            try:
                df[col] = pd.to_datetime(df[col])
            except Exception:
                pass

        if is_datetime64_any_dtype(df[col]):
            df[col] = df[col].dt.tz_localize(None)

    modification_container = st.container()

    with modification_container:
        to_filter_columns = st.multiselect("Filter dataframe on", df.columns)
        for column in to_filter_columns:
            _, right = st.columns((1, 20))
            # Treat columns with < 10 unique values as categorical
            if is_categorical_dtype(df[column]) or df[column].nunique() < 3:
                user_cat_input = right.multiselect(
                    f"Values for {column}",
                    df[column].unique(),
                    default=list(df[column].unique()),
                )
                df = df[df[column].isin(user_cat_input)]
            elif is_numeric_dtype(df[column]):
                if all([isinstance(item, (int, np.uint)) for item in df[column]]):
                    _min = int(df[column].min())
                    _max = int(df[column].max())
                    step = 1
                else:
                    _min = float(df[column].min())
                    _max = float(df[column].max())
                    step = (_max - _min) / 100

                user_num_input = right.slider(
                    f"Values for {column}",
                    min_value=_min,
                    max_value=_max,
                    value=(_min, _max),
                    step=step,)
                df = df[df[column].between(*user_num_input)]
            elif is_datetime64_any_dtype(df[column]):
                user_date_input = right.date_input(
                    f"Values for {column}",
                    value=(
                        df[column].min(),
                        df[column].max(),
                    ),
                )
                if len(user_date_input) == 2:
                    user_date_input = tuple(map(pd.to_datetime, user_date_input))
                    start_date, end_date = user_date_input
                    df = df.loc[df[column].between(start_date, end_date)]
            else:
                user_text_input = right.text_input(
                    f"Substring or regex in {column}",
                )
                if user_text_input:
                    df = df[df[column].astype(str).str.contains(user_text_input)]

    return df


# TODO use the selection of the table and download The docking pose storage in the Individual
# Or something like make_sdf of moldrug and download the info
# Another tab that print the general information how went the rum print some convergency
# The violin plot


def convert(number):
    if isinstance(number, np.floating):
        return float(number)
    elif isinstance(number, np.integer):
        return int(number)
    # It is not possible to convert, this fix some issue on MacOS
    else:
        return number


def plot_dist(individuals: list[utils.Individual], properties: list[str], every_gen: int = 1, figsize=(25, 25)):
    """Create the violin plot for the MolDrug run

    Parameters
    ----------
    individuals : list[utils.Individual]
        A list of individuals
    properties : list[str]
        A list of the properties to be graph (must be attributes of the provided individuals)
    every_gen : int, optional
        Frequency to plot the distribution: every how many generations, by default 1
    fig_size : tuple(int), optional
        The size of the graph, by default (25, 25)

    Returns
    -------
    tuple
        fig, axes
    """
    # Set up the matplotlib figure
    if len(properties) <= 1:
        extra_plot_kwargs = dict(sharex=True, gridspec_kw={'hspace': 0.05})
    else:
        extra_plot_kwargs = dict()
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(nrows=len(properties), figsize=figsize, **extra_plot_kwargs)

    SawIndividuals = utils.to_dataframe(individuals).drop(['pdbqt'], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    SawIndividuals = SawIndividuals[SawIndividuals['kept_gens'].map(len) != 0].reset_index(drop=True)
    gen_idxs = sorted(set(item for sublist in SawIndividuals['kept_gens'] for item in sublist))
    NumGens = max(gen_idxs)

    # Set pop to the initial population and pops out the first gen
    pop = SawIndividuals[SawIndividuals.genID == gen_idxs.pop(0)].sort_values(by=["cost"])
    pops = pop.copy()
    for gen_idx in gen_idxs:
        idx = [i for i in range(SawIndividuals.shape[0]) if gen_idx in SawIndividuals.loc[i, 'kept_gens']]
        pop = SawIndividuals.copy().iloc[idx, :].assign(genID=gen_idx)
        pops = pd.concat([pops, pop.copy()])
    # Draw a violinplot with a narrow bandwidth than the default
    pops = pops.loc[pops['genID'].isin([gen for gen in range(0, NumGens + every_gen, every_gen)])]

    for i, prop in enumerate(properties):
        if len(properties) <= 1:
            sns.violinplot(hue='genID', x='genID', y=prop, data=pops, palette="Set3", bw_adjust=.2, cut=0, linewidth=1, ax=axes, legend=False)
        else:
            sns.violinplot(hue='genID', x='genID', y=prop, data=pops, palette="Set3", bw_adjust=.2, cut=0, linewidth=1, ax=axes[i], legend=False)

            # Remove x-axis labels and ticks from all but the bottom subplot
            if i < len(properties) - 1:
                axes[i].xaxis.set_visible(False)
    return fig, axes



@st.cache_data
def ProtPdbBlockToProlifMol(protein_pdb_string):
    with tempfile.NamedTemporaryFile(prefix='.pro', suffix='.pdb', mode='w+') as tmp:
        tmp.write(protein_pdb_string)
        protein = mda.Universe(tmp.name)
        protein = plf.Molecule.from_mda(protein)
    return protein


def MolFromPdbqtBlock(pdbqt_string):
    pdbqt_tmp = tempfile.NamedTemporaryFile(suffix='.pdbqt')
    with open(pdbqt_tmp.name, 'w') as f:
        f.write(pdbqt_string)
    pdbqt_mol = PDBQTMolecule.from_file(pdbqt_tmp.name, skip_typing=True)
    mol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]
    return mol


def LigPdbqtBlockToProlifMol(ligand_pdbqt_string):
    ligand = MolFromPdbqtBlock(ligand_pdbqt_string)
    ligand = plf.Molecule.from_rdkit(ligand)
    return ligand


@st.cache_data
def prolif_plot_2d(ligand_pdbqt_string, protein_pdb_string):

    # ProLIF example
    # load topology
    # Protein
    protein = ProtPdbBlockToProlifMol(protein_pdb_string)
    ligand = LigPdbqtBlockToProlifMol(ligand_pdbqt_string)

    fp = plf.Fingerprint()
    fp.run_from_iterable([ligand], protein)
    try:
        prolif_ligplot_html_document = fp.plot_lignetwork(
            ligand,
            # replace with `kind="frame", frame=0` for the other depiction
            kind="aggregate",
            threshold=0.3,
            rotation=270,
            height="400px"
        ).data

    except ValueError:
        prolif_ligplot_html_document = None
    return prolif_ligplot_html_document


@st.cache_data
def prolif_plot_3d(ligand_pdbqt_string, protein_pdb_string, spin=False):

    # ligand = MolFromPdbqtBlock(ligand_pdbqt_string)
    # view = py3Dmol.view()
    # view.removeAllModels()
    # view.setViewStyle({'style':'outline','color':'black','width':0.1})

    # view.addModel(protein_pdb_string,format='pdb')
    # Prot=view.getModel()
    # Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})
    # view.addSurface(py3Dmol.VDW,{'opacity':0.6,'color':'white'})

    # view.addModel(Chem.MolToMolBlock(ligand),format='mol2')
    # ref_m = view.getModel()
    # ref_m.setStyle({},{'stick':{'colorscheme':'greenCarbon','radius':0.2}})
    # if spin:
    #     view.spin(True)
    # else:
    #     view.spin(False)
    # view.zoomTo()

    protein = ProtPdbBlockToProlifMol(protein_pdb_string)
    ligand = LigPdbqtBlockToProlifMol(ligand_pdbqt_string)

    fp = plf.Fingerprint()
    fp.run_from_iterable([ligand], protein)
    try:
        view = fp.plot_3d(
            ligand_mol=ligand,
            protein_mol=protein,
            frame=0,
            # replace with `kind="frame", frame=0` for the other depiction
            display_all=False,
        )
        if spin:
            view.spin(True)
        else:
            view.spin(False)
        showmol(view, height=500, width=800)
    except ValueError:
        st.warning('No possible to display')


@st.cache_data
def lig_prot_overview(_pop, protein_pdb_string):
    protein = ProtPdbBlockToProlifMol(protein_pdb_string)
    with tempfile.NamedTemporaryFile(mode='w+') as tmp:
        utils.make_sdf(_pop, tmp.name)
        directory, name = os.path.split(tmp.name)
        lig_suppl = plf.sdf_supplier(os.path.join(directory, f"{name}.sdf"))
    # generate fingerprint
    fp = plf.Fingerprint()
    fp.run_from_iterable(lig_suppl, protein)
    df = fp.to_dataframe()
    df = df.droplevel("ligand", axis=1)

    # aminoacids = set()
    # interaction_types = set()
    # for aminoacid, interaction_type in df.columns:
    #     aminoacids.add(aminoacid)
    #     interaction_types.add(interaction_type)
    # aminoacids = sorted(aminoacids)
    # interaction_types = sorted(interaction_types)

    df.index = [individual.idx for individual in _pop]
    df.index.names = ['idx']
    df = df.groupby(level='interaction', axis=1).sum()

    # show all interactions with a specific protein residue
    # df = df.xs("TYR31.A", level="protein", axis=1)
    # show all pi-stacking interactions
    # df = df.xs(interaction_types[0], level="interaction", axis=1)
    return df


@st.cache_data
def load_pbz2(pbz2):
    moldrug_result = utils.decompress_pickle(pbz2)
    is_GA = False
    if isinstance(moldrug_result, utils.GA):
        gen, pop = moldrug_result.NumGens, moldrug_result.pop
        is_GA = True
    elif isinstance(moldrug_result, utils.Local):
        gen, pop = 0, moldrug_result.pop
    elif isinstance(moldrug_result, tuple):
        if isinstance(moldrug_result[0], int) and isinstance(moldrug_result[1][0], utils.Individual):
            gen, pop = moldrug_result[0], moldrug_result[1]
    else:
        raise Exception('pbz2 is corrupted')

    try:
        dataframe = utils.to_dataframe(pop, return_mol=True)
    except TypeError:
        dataframe = utils.to_dataframe(pop)
    pdbqt_dataframe = dataframe[['idx', 'pdbqt']]
    pdbqt_dataframe.set_index('idx', inplace=True)
    return gen, moldrug_result, pdbqt_dataframe, dataframe, is_GA


@st.cache_data
def upload_file_to_string(uploaded_file):
    # To convert to a string based IO:
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    # To read file as string:
    string_data = stringio.read()
    return string_data


# PubChem functions and interaction with web servers
def get_chemazone_price(smiles):
    encoded_smiles = urllib.parse.quote(smiles)

    url = f"https://chemazone.com/StrSearch.asp?allfields={encoded_smiles}"
    # Make a POST request with the SMILES data
    response = requests.get(url)

    # Parse the HTML content using BeautifulSoup
    soup = BeautifulSoup(response.content, "html.parser")
    # Find the relevant HTML element
    option_tag = soup.find("option")
    try:
        return option_tag['value'], option_tag['data-price']
    except TypeError:
        return '-', '-'


def get_similarity(smiles1: str, smiles2: str) -> float:
    """Get the similarity between two molecules

    Parameters
    ----------
    smiles1 : str
        The SMILES representation of the first molecule
    smiles2 : str
        The SMILES representation of the second molecule

    Returns
    -------
    float
        Similarity
    """
    # Convert SMILES strings to RDKit molecules
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    # Calculate fingerprints for each molecule

    fp1 = Chem.RDKFingerprint(mol1)
    fp2 = Chem.RDKFingerprint(mol2)

    # Calculate similarity
    similarity = DataStructs.FingerprintSimilarity(fp1, fp2)

    return similarity


@st.cache_data
def get_compound_vendors(cid: int) -> tuple:
    """Gte the vendors of the molecule with CID in PubChem

    Parameters
    ----------
    cid : int
        Compound Identification in PubChem

    Returns
    -------
    tuple
        (vendor URL; number of vendors)
    """
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/{cid}/JSON'

    response = requests.get(url)

    data = response.json()
    if 'SourceCategories' in data:
        num_sources = len(data['SourceCategories']['Categories'])
        vendor_url = data['SourceCategories']['Categories'][0]['URL']
        return vendor_url, num_sources
    else:
        return None, 0


@st.cache_data
def get_pubchem_data(smiles: str, Threshold: int = 95) -> dict:
    """Get the data from PubChem

    Parameters
    ----------
    smiles : str
        The SMILES identification
    Threshold : int, optional
        Threshold of similarity in case that the molecule is not in PubChem, by default 95

    Returns
    -------
    dict
        A dictionary with keywords: smiles, cid, similarity, vendors_link, num_vendors
    """
    result = {
        'smiles': None,
        'cid': None,
        'similarity': 0,
        'vendors_link': None,
        'num_vendors': 0,
        # 'chemazone_price': None,
    }
    num_iter = 1
    onsimilarity = False
    for i in range(num_iter):
        try:
            if not onsimilarity:
                compound = pcp.get_compounds(
                    identifier=smiles,
                    namespace='smiles',
                    domain='compound')[0]
                if compound.cid:
                    vendors_link, num_vendors = get_compound_vendors(compound.cid)
                    result.update({
                        'smiles': smiles,
                        'cid': compound.cid,
                        'similarity': 1,
                        'vendors_link': vendors_link,
                        'num_vendors': num_vendors,
                        # 'chemazone_price': '/'.join(get_chemazone_price(smiles))
                    })
                    break
                else:
                    onsimilarity = True
            if onsimilarity:
                compound = pcp.get_compounds(
                    identifier=smiles,
                    namespace='smiles',
                    domain='compound',
                    searchtype='similarity',
                    Threshold=Threshold,
                    MaxRecords=1)[0]
                if compound.cid:
                    vendors_link, num_vendors = get_compound_vendors(compound.cid)
                    result.update({
                        'smiles': compound.isomeric_smiles,
                        'cid': compound.cid,
                        'similarity': get_similarity(smiles, compound.isomeric_smiles),
                        'vendors_link': vendors_link,
                        'num_vendors': num_vendors,
                        # 'chemazone_price': '/'.join(get_chemazone_price(compound.isomeric_smiles)),
                    })
                    break
        except pcp.PubChemHTTPError as e:
            if e.msg == 'PUGREST.ServerBusy':
                time.sleep(3)
                if i == num_iter - 1:
                    result.update({'cid': e.msg})
            else:
                result.update({'cid': e.msg})
                break

    return result


@st.cache_data
def get_pubchem_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """It calls get_pubchem_data on each smiles of df and return
    a DataFrame

    Parameters
    ----------
    df : pd.DataFrame
        With columns: idx (the identification used in MolDrug) and smiles (the SMILES of the molecule)

    Returns
    -------
    pd.DataFrame
        _description_
    """
    data = []
    for _, row in df.iterrows():
        new_row = get_pubchem_data(smiles=row['smiles'])
        new_row['idx'] = row['idx']
        data.append(new_row)
    return pd.DataFrame(data)[['idx', 'cid', 'similarity', 'smiles', 'num_vendors', 'vendors_link']]  # , 'chemazone_price']]

    # data = []
    # with concurrent.futures.ThreadPoolExecutor() as executor:
    #     future_to_row = {executor.submit(get_pubchem_data, row['smiles']): row for _, row in df.iterrows()}
    #     for future in concurrent.futures.as_completed(future_to_row):
    #         row = future_to_row[future]
    #         new_row = future.result()
    #         new_row['idx'] = row['idx']
    #         data.append(new_row)
    # return pd.DataFrame(data)[['idx', 'cid', 'similarity', 'smiles', 'num_vendors', 'vendors_link']]


if __name__ == "__main__":
    # Upload the data, result and PDB used fopr the docking
    st.sidebar.subheader('Upload pbz2:')
    # pbz2 = '/home/ale/mnt/snowden2/MolDrug/HIPS/Second/LasB/Cost/free/CH4/replica1/04_local_result.pbz2'
    pbz2 = st.sidebar.file_uploader('**MolDrug pbz2**', accept_multiple_files=False)

    if pbz2:
        # TODO, make gen a selectable user parameter
        gen, moldrug_result, pdbqt_dataframe, dataframe, is_GA = load_pbz2(pbz2)
        st.sidebar.write(f"**Loaded generation = {gen}**")
        try:
            dataframe['mol'] = dataframe['mol'].apply(Chem.RemoveHs)
        except KeyError:
            dataframe['mol'] = dataframe['pdbqt'].apply(lambda x: Chem.RemoveHs(MolFromPdbqtBlock(x)))
        properties = [prop for prop in dataframe.columns if prop not in ['idx', 'pdbqt', 'mol', 'kept_gens']]
        properties = st.sidebar.multiselect("Choose properties", properties, ["cost"])
        # Plot the prolif and the tridimensional structure with py3Dmol
        # That two columns
        # Get the minimum and maximum of the variables and used them
        min_cost, max_cost = convert(dataframe['cost'].min()), convert(dataframe['cost'].max())
        cost_threshold = st.sidebar.slider('**Coloring cost threshold:**', min_cost, max_cost + 0.001, (min_cost+max_cost)/2)
        sliders = []
        for prop in properties:
            minimum, maximum = convert(dataframe[prop].min()), convert(dataframe[prop].max())
            if minimum == maximum:
                maximum += 0.001
            sliders.append(st.sidebar.slider(f'**{prop}**', minimum, maximum, [minimum, maximum]))

        grid = mols2grid.MolGrid(dataframe, mol_col='mol', fixedBondLength=25, clearBackground=False, size=(130, 120))

        for prop, slide in zip(properties, sliders):
            grid.dataframe = grid.dataframe[(slide[0] <= grid.dataframe[prop]) & (grid.dataframe[prop] <= slide[1])]

        try:
            view = grid.display(
                # set what's displayed on the grid
                subset=["idx", "img", "cost"],
                # # set what's displayed on the hover tooltip
                tooltip=list(set(["cost", 'idx'] + properties)),
                tooltip_placement="auto",
                # style for the grid labels and tooltips
                style={
                    "cost": lambda x: "color: red; font-weight: bold;" if x < cost_threshold else "",
                    "__all__": lambda x: "background-color: azure;" if x["cost"] >= cost_threshold else ""
                },
                # change the precision and format (or other transformations)
                transform={"cost": lambda x: round(x, 3)},
                # sort the grid in a different order by default
                sort_by="cost",
                n_rows=3,
                callback=mols2grid.callbacks.info(img_size=(200, 150)),
            )

            with tab1:
                components.html(view.data, width=None, height=700, scrolling=True)
                with st.expander('Show table of properties'):
                    # For compatibility with older versions
                    props_to_drop = ['mol', 'pdbqt', 'img', 'mols2grid-id']
                    if 'kept_gens' in grid.dataframe.columns:
                        props_to_drop.append('kept_gens')

                    prop_df = grid.dataframe.drop(['mol', 'pdbqt', 'kept_gens', 'img', 'mols2grid-id'], axis=1).set_index('idx')
                    st.download_button(
                        "Press to Download",
                        convert_df(prop_df),
                        "MolDrug_properties.csv",
                        "text/csv",
                        key='download-MolDrug_prop-csv')
                    st.dataframe(prop_df)

        except ValueError:
            with tab1:
                st.info('Nothing to show')

        with tab3:
            st.info('🔬 Please note that this feature is currently in the experimental phase. '
                    'Occasionally, queries may encounter issues and fail to generate the expected results. 🚧')
            PubChemCheck = st.checkbox("Explore PubChem")
            # TODO Use https://github.com/whitead/molbloom as well
            download_button = st.empty()
            if PubChemCheck:
                dataframe = dataframe.copy()
                dataframe['smiles'] = dataframe['mol'].apply(Chem.MolToSmiles)
                # It consumes to much resources, for the moment, we will include some manual index selections
                # dataframe['idx'][:2]
                idx_selection = st.multiselect(label="Choose some idxs", options=grid.dataframe['idx'], default=None)
                if idx_selection:
                    pubchem_dataframe = get_pubchem_dataframe(dataframe[dataframe[['idx', 'smiles']]['idx'].isin(idx_selection)])
                    # pubchem_dataframe = pubchem_dataframe[pubchem_dataframe['idx'].isin(grid.dataframe['idx'])]
                    pubchem_dataframe = pubchem_dataframe.set_index('idx')
                    st.dataframe(pubchem_dataframe)
                    download_button.download_button(
                        "Press to Download",
                        convert_df(pubchem_dataframe),
                        "PubChemData.csv",
                        "text/csv",
                        key='download-PubChemData-csv')
                else:
                    st.info('☝️ You must select something. Be patient, this could take a while ⌛')

        plif = st.empty()

        st.sidebar.subheader('**Ligand-protein network interaction**')

        # Every form must have a submit button
        protein_pdb = st.sidebar.file_uploader('**Protein PDB**', accept_multiple_files=False)
        if protein_pdb:
            protein_pdb_string = upload_file_to_string(protein_pdb)

            # Get overview
            if is_GA:
                df_overview = lig_prot_overview(moldrug_result.pop, protein_pdb_string=protein_pdb_string)
            else:
                df_overview = lig_prot_overview(moldrug_result[1], protein_pdb_string=protein_pdb_string)

            # Input widget
            idx = st.sidebar.selectbox('idx', sorted(grid.dataframe['idx']))
            representation = st.sidebar.selectbox('Representation', ['2D', '3D'])
            spin = st.sidebar.empty()
            if representation == '2D':
                prolif_ligplot_html = prolif_plot_2d(
                    ligand_pdbqt_string=pdbqt_dataframe.loc[idx, 'pdbqt'],
                    protein_pdb_string=protein_pdb_string,
                )
                with plif:
                    with tab1:
                        if prolif_ligplot_html:
                            components.html(prolif_ligplot_html, width=None, height=500, scrolling=True)
                        else:
                            st.error("It was not possible to genereate the ProLIF image.")
            else:
                spin = spin.checkbox('Spin', value=False)
                with plif:
                    with tab1:
                        prolif_plot_3d(
                            ligand_pdbqt_string=pdbqt_dataframe.loc[idx, 'pdbqt'],
                            protein_pdb_string=protein_pdb_string,
                            spin=spin)
            with tab1:
                with st.expander('Table of interactions'):
                    filter_df_overview = filter_dataframe(df_overview[df_overview.index.isin(grid.dataframe['idx'])])
                    st.download_button(
                        "Press to Download",
                        convert_df(filter_df_overview),
                        "Protein-ligand_interactions.csv",
                        "text/csv",
                        key='download-pli-csv')
                    st.dataframe(filter_df_overview)
        else:
            st.sidebar.info('☝️ Upload the PDB protein file.')

        # Plot the distribution
        with tab2:
            try:
                every_gen = st.number_input("Every how many generations:", min_value=1, max_value=moldrug_result.NumGens, value=10)
                col1, col2 = st.columns(2)

                # Add widgets to each column
                fig_size_x = col1.number_input("Size of the figure in x:", value=10)
                fig_size_y = col2.number_input("Size of the figure in y:", value=10)
                properties_to_plot = [prop for prop in properties if prop not in ['genID']]
                fig, axes = plot_dist(moldrug_result.SawIndividuals, properties=properties_to_plot, every_gen=every_gen, figsize=(fig_size_x, fig_size_y))
                st.pyplot(fig)
            except Exception:
                if is_GA:
                    st.info('Nothing to show. Consider to select some properties in the side bar.')
                else:
                    st.info('Nothing to show. The input is not a MoDrug GA.')

    else:
        st.sidebar.info("☝️ Upload MolDrug's pbz2 file.")
