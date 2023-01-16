import mols2grid
from rdkit import Chem
from prolif.plotting.network import LigNetwork
import prolif as plf
import MDAnalysis as mda
import numpy as np
import streamlit as st
import tempfile
import streamlit.components.v1 as components
from moldrug import utils
from meeko import RDKitMolCreate, PDBQTMolecule
from io import StringIO
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import py3Dmol
from stmol import showmol

# st.set_page_config('wide')
st.sidebar.empty()
st.title('Dashboard')
st.image('https://github.com/ale94mleon/MolDrug/raw/main/docs/source/_static/logo.png?raw=true', width=150)

with st.expander('**About the App**'):
    st.markdown("👈 Open the side bar to introduce the data.\n\n"\
        "This app is to get an overview of a Moldrug result at glance. \n"
        "Check [MolDrug's docs](https://moldrug.rtfd.io/) and [MolDrug's GitHub](https://github.com/ale94mleon/moldrug/) for more information.")

tab1, tab2 = st.tabs(["Molecules", "Running info"])



def MolFromPdbqtBlock(pdbqt_string):
    pdbqt_tmp = tempfile.NamedTemporaryFile(suffix='.pdbqt')
    with open(pdbqt_tmp.name, 'w') as f:
        f.write(pdbqt_string)
    pdbqt_mol = PDBQTMolecule.from_file(pdbqt_tmp.name, skip_typing=True)
    mol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]
    return mol

#TODO use the selection of the table and download The docking pose storage in the Individual
# Or something like make_sdf of moldrug and download the info
# Another tab that print the general information how went the rum print some convergency
# The violin plot

def convert(number):
    if isinstance(number, np.floating):
        return float(number)
    if isinstance(number, np.integer):
        return int(number)

def plot_dist(individuals:list[utils.Individual], properties:list[str], every_gen:int = 1):
    """Create the violin plot for the MolDrug run

    Parameters
    ----------
    individuals : list[utils.Individual]
        A list of individuals
    properties : list[str]
        A list of the properties to be graph (must be attributes of the provided individuals)
    every_gen : int, optional
        Frequency to plot the distrubution: every how many generations, by default 1

    Returns
    -------
    tuple
        fig, axes
    """


    # Set up the matplotlib figure
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(nrows = len(properties), figsize=(25, 25))

    SawIndividuals = utils.to_dataframe(individuals).drop(['pdbqt'], axis = 1).replace([np.inf, -np.inf], np.nan).dropna()
    SawIndividuals = SawIndividuals[SawIndividuals['kept_gens'].map(len) != 0].reset_index(drop=True)
    gen_idxs = sorted(SawIndividuals.genID.unique())
    NumGens = max(gen_idxs)

    # Set pop to the initial population and pops out the first gen
    pop = SawIndividuals[SawIndividuals.genID == gen_idxs.pop(0)].sort_values(by=["cost"])
    pops = pop.copy()
    for gen_idx in gen_idxs:
        idx = [i for i in range(SawIndividuals.shape[0]) if gen_idx in SawIndividuals.loc[i,'kept_gens']]
        pop = SawIndividuals.copy().iloc[idx,:].assign(genID=gen_idx)
        pops = pd.concat([pops, pop.copy()])
    # Draw a violinplot with a narrow bandwidth than the default
    pops = pops.loc[pops['genID'].isin([gen for gen in range(0, NumGens+every_gen, every_gen)])]

    if len(properties) <= 1:
        sns.violinplot(x = 'genID', y = properties[0], data=pops, palette="Set3", bw=.2, cut=0, linewidth=1, ax=axes)
    else:
        for i, prop in enumerate(properties):
            sns.violinplot(x = 'genID', y = prop, data=pops, palette="Set3", bw=.2, cut=0, linewidth=1, ax=axes[i])

    return fig, axes

@st.cache
def ProtPdbBlockToProlifMol(protein_pdb_string):
    with tempfile.NamedTemporaryFile(prefix='.pro', suffix='.pdb', mode='w+') as tmp:
        tmp.write(protein_pdb_string)
        protein = mda.Universe(tmp.name)
        protein = plf.Molecule.from_mda(protein)
    return protein

@st.cache
def LigPdbqtBlockToProlifMol(ligand_pdbqt_string):
    ligand = MolFromPdbqtBlock(ligand_pdbqt_string)
    ligand = plf.Molecule.from_rdkit(ligand)
    return ligand

@st.cache
def prolif_plot(ligand_pdbqt_string,protein_pdb_string):

    # ProLIF example
    # load topology
    # Protein
    protein = ProtPdbBlockToProlifMol(protein_pdb_string)
    ligand = LigPdbqtBlockToProlifMol(ligand_pdbqt_string)

    fp = plf.Fingerprint()
    fp.run_from_iterable([ligand], protein)
    df_fp = fp.to_dataframe(return_atoms=True)

    net = LigNetwork.from_ifp(
        df_fp,
        ligand,
        # replace with `kind="frame", frame=0` for the other depiction
        kind="aggregate",
        threshold=0.3,
        rotation=270,
    )

    prolif_ligplot_html_document = net.display(height="400px").data
    return prolif_ligplot_html_document

def py3Dmol_plot(ligand_pdbqt_string,protein_pdb_string, spin = False):

    ligand = MolFromPdbqtBlock(ligand_pdbqt_string)
    view = py3Dmol.view()
    view.removeAllModels()
    view.setViewStyle({'style':'outline','color':'black','width':0.1})

    view.addModel(protein_pdb_string,format='pdb')
    Prot=view.getModel()
    Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})
    view.addSurface(py3Dmol.VDW,{'opacity':0.6,'color':'white'})

    view.addModel(Chem.MolToMolBlock(ligand),format='mol2')
    ref_m = view.getModel()
    ref_m.setStyle({},{'stick':{'colorscheme':'greenCarbon','radius':0.2}})
    if spin:
        view.spin(True)
    else:
        view.spin(False)
    view.zoomTo()
    showmol(view,height=500,width=800)


@st.cache(allow_output_mutation=True)
def load_pbz2(pbz2):
    moldrug_result = utils.decompress_pickle(pbz2)
    if isinstance(moldrug_result, utils.GA):
        gen, pop = moldrug_result.NumGens, moldrug_result.pop
    elif isinstance(moldrug_result, utils.Local):
        gen, pop = 0, moldrug_result.pop
    elif isinstance(moldrug_result, tuple):
        if isinstance(moldrug_result[0], int) and isinstance(moldrug_result[1][0], utils.Individual):
            gen, pop = moldrug_result[0], moldrug_result[1]
    else:
        raise Exception('pbz2 is corrupted')


    try:
        dataframe = utils.to_dataframe(pop, return_mol = True)
    except TypeError:
        dataframe = utils.to_dataframe(pop)
    pdbqt_dataframe = dataframe[['idx','pdbqt']]
    pdbqt_dataframe.set_index('idx', inplace=True)
    return gen, moldrug_result, pdbqt_dataframe, dataframe

@st.cache
def upload_file_to_string(uploaded_file):
    # To convert to a string based IO:
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    # To read file as string:
    string_data = stringio.read()
    return string_data


# Upload the data, result and PDB used fopr the docking
st.sidebar.subheader('Upload pbz2:')
# pbz2 = '/home/ale/mnt/snowden2/MolDrug/HIPS/Second/LasB/Cost/free/CH4/replica1/04_local_result.pbz2'
pbz2 = st.sidebar.file_uploader('**MolDrug pbz2**', accept_multiple_files = False)



if pbz2:

    gen, moldrug_result, pdbqt_dataframe, dataframe = load_pbz2(pbz2)
    st.sidebar.write(f"**Loaded generation = {gen}**")
    try:
        dataframe['mol'] = dataframe['mol'].apply(Chem.RemoveHs)
    except KeyError:
        dataframe['mol'] = dataframe['pdbqt'].apply(lambda x: Chem.RemoveHs(MolFromPdbqtBlock(x)))

    properties = [prop for prop in dataframe.columns if prop not in ['idx','pdbqt', 'mol', 'kept_gens']]
    properties = st.sidebar.multiselect(
    "Choose properties", properties, ["cost"]
    )


    # # Plot the prolif and the tridimensional structure with py3Dmol
    # # That two columns


    # Get the minimum and maximum of the variables and used them
    cost_threshold = st.sidebar.slider('**Coloring cost threshold:**', 0.0,1.0,0.5)
    sliders = []
    for prop in properties:
        minimum, maximum = convert(dataframe[prop].min()), convert(dataframe[prop].max())
        sliders.append(st.sidebar.slider(f'**{prop}**', minimum,maximum,[minimum,maximum]))


    grid = mols2grid.MolGrid(
        dataframe,
        mol_col = 'mol',
        # molecule drawing parameters
        fixedBondLength=25,
        clearBackground=False,
        size=(130, 120),
    )

    for prop, slide in zip(properties,sliders):
        grid.dataframe = grid.dataframe[(slide[0]<= grid.dataframe[prop]) & (grid.dataframe[prop]<=slide[1])]

    try:
        view = grid.display(
            # set what's displayed on the grid
            subset=["idx", "img", "cost"],
            # # set what's displayed on the hover tooltip
            tooltip=list(set(["cost",'idx'] + properties)),
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

    except ValueError:
        with tab1:
            st.info('Nothing to show')


    plif = st.empty()

    st.sidebar.subheader('**Ligand-protein network interaction**')


    # Every form must have a submit button
    protein_pdb = st.sidebar.file_uploader('**Protein PDB**', accept_multiple_files = False)
    if protein_pdb:
        protein_pdb_string=upload_file_to_string(protein_pdb)
        # Input widge
        idx = st.sidebar.selectbox('idx', sorted(pdbqt_dataframe.index))
        representation = st.sidebar.selectbox('Representation', ['ProLIF', 'Py3Dmol'])
        spin = st.sidebar.empty()
        if representation == 'ProLIF':
            prolif_ligplot_html_document = prolif_plot(
                ligand_pdbqt_string=pdbqt_dataframe.loc[idx, 'pdbqt'],
                protein_pdb_string=protein_pdb_string,
            )
            with plif:
                with tab1:
                    components.html(prolif_ligplot_html_document,width=None, height=500, scrolling=True)
        else:
            spin.checkbox('Spin', value = False)
            with plif:
                with tab1:
                    py3Dmol_plot(
                        ligand_pdbqt_string=pdbqt_dataframe.loc[idx, 'pdbqt'],
                        protein_pdb_string=protein_pdb_string,
                        spin = spin
                    )
    else:
        st.sidebar.info('☝️ Upload the PDB protein file.')

    # Plot the distribution
    with tab2:
        every_gen = st.number_input("Every how many generations:",min_value=1, max_value=moldrug_result.NumGens, value=10)
        try:
            properties_to_plot = [prop for prop in properties if prop not in ['genID']]
            fig, axes = plot_dist(moldrug_result.SawIndividuals,properties=properties_to_plot, every_gen=every_gen)
            st.pyplot(fig)
        except Exception as e:
            st.info('Nothing to show. Consider to select some properties in the side bar.')

else:
    st.sidebar.info('☝️ Upload a MolDrug pbz2 file.')