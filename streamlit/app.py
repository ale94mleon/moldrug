import mols2grid
from rdkit import Chem
from prolif.plotting.network import LigNetwork
import prolif as plf
import numpy as np
import streamlit as st
import tempfile
import streamlit.components.v1 as components
from moldrug import utils
from meeko import RDKitMolCreate, PDBQTMolecule
from io import StringIO

# st.set_page_config('wide')
st.sidebar.empty()
st.title('Dashboard')
st.image('https://github.com/ale94mleon/MolDrug/raw/main/docs/source/_static/logo.png?raw=true', width=150)

with st.expander('**More Info...**'):
    st.markdown("Check the [MolDrug](https://moldrug.rtfd.io/) documentation or our [GitHub Repo](https://github.com/ale94mleon/moldrug/).")


def MolFromPdbqtBlock(pdbqt_string):
    pdbqt_tmp = tempfile.NamedTemporaryFile(suffix='.pdbqt')
    with open(pdbqt_tmp.name, 'w') as f:
        f.write(pdbqt_string)
    pdbqt_mol = PDBQTMolecule.from_file(pdbqt_tmp.name, skip_typing=True)
    mol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]
    return mol

#TODO use the selction of the table and download The docking pose storge in the Individual
# Or something like make_sdf of moldrug and download the info
# Another tab that print the general information how went the rum print some convergency
# The violin plot

def convert(number):
    if isinstance(number, np.floating):
        return float(number)
    if isinstance(number, np.integer):
        return int(number)


def prolif_plot(ligand_pdbqt_string,protein_pdb_string):

    # ProLIF example
    # load topology
    # Protein
    protein = Chem.MolFromPDBBlock(protein_pdb_string, removeHs = False)
    protein = plf.Molecule.from_rdkit(protein)
    ligand = MolFromPdbqtBlock(ligand_pdbqt_string)
    ligand = plf.Molecule.from_rdkit(ligand)
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

    # prolif_plot_callback_return = mols2grid.make_popup_callback(
    #     title="${data['SMILES']}",
    #     html=prolif_ligplot_html_document,
    #     style="max-width: 80%;"
    # )
    return prolif_ligplot_html_document



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
    with st.spinner("Recalculating..."):

        gen, _, pdbqt_dataframe, dataframe = load_pbz2(pbz2)
        st.sidebar.write(f"**Loaded generation = {gen}**")
        try:
            dataframe['mol'] = dataframe['mol'].apply(lambda x: Chem.RemoveHs(x))
        except KeyError:
            dataframe['mol'] = dataframe['pdbqt'].apply(lambda x: Chem.RemoveHs(MolFromPdbqtBlock(x)))

        properties = [prop for prop in dataframe.columns if prop not in ['idx','pdbqt', 'mol', 'kept_gens']]
        properties = st.multiselect(
        "Choose properties", properties, ["cost"]
        )
        # # Identify the variables
        # # Create a sidebar with the variables from the minimum to the maximum, or the user must
        # # introduce the name of the variables used, more or less select from the possible variables
        # # all of this in the sidebar
        # # Plot the mols2grid, I copuld even try with a normal dataframe and see if is acceptable
        # # Plot the prolif and the tridementionl structure with py3Dmol
        # # That two columns


        # Get the minimum and maximum of the variables and used them
        cost_threshold = st.sidebar.slider('**Coloring cost threshold:**', 0.0,1.0,0.5)
        sliders = []
        for prop in properties:
            minimum, maximum = convert(dataframe[prop].min()), convert(dataframe[prop].max())
            sliders.append(st.sidebar.slider(f'**{prop}**', minimum,maximum,[minimum,maximum]))        


        # for i, prop in enumerate(properties):
        #     minimum, maximum = dataframe[prop].min()+1, dataframe[prop].max()
        #     st.write(prop, minimum, maximum)
        # cost = st.sidebar.slider('**cost**', 0.0,1.0,[0.0,1.0])
        # vina_score = st.sidebar.slider('**vina_score**', -13.0,5.0,[-13.0,5.0])
        # sa_score = st.sidebar.slider('**sa_score**', 0.0,13.0,[0.0,13.0])
        # qed = st.sidebar.slider('**qed**', 0.0,1.0,[0.0,1.0])
        # query = "@cost[0] <= cost <= @cost[1] and "\
        #     "@vina_score[0] <= vina_score <= @vina_score[1] and "\
        #     "@sa_score[0] <= sa_score <= @sa_score[1] and "\
        #     "@qed[0] <= qed <= @qed[1]"\

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
            row1_1, row1_2 = st.columns((2,3))
            selection = mols2grid.get_selection()
            
            components.html(view.data, width=None, height=700, scrolling=True)
            row1_1, row1_2 = st.columns((2,3))
        
        except ValueError:
            st.info('Nothing to show')
        
        
        plif = st.empty()
        
        st.sidebar.subheader('**Ligand-protein network interacion**')


        # Every form must have a submit button
        protein_pdb = st.sidebar.file_uploader('**Protein PDB**', accept_multiple_files = False)
        if protein_pdb:
            # Input widge
            idx = st.sidebar.selectbox('idx', sorted(pdbqt_dataframe.index))
            representation = st.sidebar.selectbox('Representation', ['ProLIF', 'Py3Dmol'])
            if representation == 'ProLIF':
                protein_pdb_string=upload_file_to_string(protein_pdb)
                
                prolif_ligplot_html_document = prolif_plot(
                    ligand_pdbqt_string=pdbqt_dataframe.loc[idx, 'pdbqt'],
                    protein_pdb_string=protein_pdb_string,
                )
                with plif:
                    components.html(prolif_ligplot_html_document,width=None, height=500, scrolling=True)
            else:
                    st.sidebar.info('Ups, it is not implemented yet')
        else:
            st.sidebar.info('☝️ Upload the PDB protein file.')

else:
    st.sidebar.info('☝️ Upload a MolDrug pbz2 file.')