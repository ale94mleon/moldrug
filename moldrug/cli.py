#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For information of MolDrug:
    Docs: https://moldrug.readthedocs.io/en/latest/
    Source Code: https://github.com/ale94mleon/moldrug
"""
from moldrug import utils, constraintconf, __version__
import yaml, argparse, inspect, os, sys, datetime
from rdkit import Chem
def __moldrug_cmd():
    """
    This function is only used in as part of the command line interface of MolDrug.
    It makes possible to use MolDrug form the command line. More detail help is available
    from the command line `moldrug -h`.

    Raises
    ------
    NotImplementedError
        In case that the type of the calculation differs from Local or GA (currently implementations)
    ValueError
        In case that the user ask for followed jobs and Local is selected.
    ValueError
        In case that a non-mutable or non-defined argument is given by the user on the follow jobs.
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        help='The configuration yaml file',
        dest='yaml_file',
        type=str)
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f"moldrug: {__version__}")
    parser.add_argument('-f', '--fitness',
                        help="The path to the user-custom fitness module; inside of which the given custom cost function must be implemented. "\
                            "See the docs for how to do it properly. E.g. my/awesome/fitness_module.py."\
                            "By default will look in the moldrug.fitness module.",
                        dest='fitness',
                        default=None,
                        type=str)
    parser.add_argument('-o', '--outdir',
                        help="The path to where all the files should be written. "\
                            "By default the current working directory will be used (where the command line was invoked).",
                        dest='outdir',
                        default=None,
                        type=str)
    args = parser.parse_args()

    with open(args.yaml_file, 'r') as c:
        Config = yaml.safe_load(c)

    # Checking how many jobs are defined: Main and follows jobs
    MainConfig = Config.pop(list(Config.keys())[0])
    FollowConfig = Config

    if args.fitness:
        # If the fitness module provided is not in the current directory or if its name is not fitness
        # Create the module inside MolDrug
        if args.outdir:
            if not os.path.exists(args.outdir): os.makedirs(args.outdir)
            destination_path = os.path.join(args.outdir, 'CustomMolDrugFitness.py')
        else:
            destination_path = 'CustomMolDrugFitness.py'
        with open(args.fitness, 'r') as source:
            with open(destination_path, 'w') as destination:
                destination.write(source.read())
        # Changing to the outdir path if provided
        if args.outdir: os.chdir(args.outdir)
        sys.path.append('.')
        import CustomMolDrugFitness
        Cost = dict(inspect.getmembers(CustomMolDrugFitness))[MainConfig['costfunc']]
    else:
        from moldrug import fitness
        Cost = dict(inspect.getmembers(fitness))[MainConfig['costfunc']]

    if MainConfig['type'].lower() == 'ga':
        TypeOfRun = utils.GA
    elif MainConfig['type'].lower() == 'local':
        TypeOfRun = utils.Local
    else:
        raise NotImplementedError(f"\"{MainConfig['type']}\" it is not a possible type. Select from: GA or Local")

    # Convert the SMILES (or path to compress_pickle) to RDKit mol (or list of RDkit mol)
    if MainConfig['type'].lower() == 'local':
        MainConfig['seed_mol'] = Chem.MolFromSmiles(MainConfig['seed_mol'])
    else:
        if isinstance(MainConfig['seed_mol'], list):
            # If the items are path to the pickle objects
            if any([os.path.isfile(path) for path in MainConfig['seed_mol']]):
                seed_pop = set()
                for solution in MainConfig['seed_mol']:
                    _, pop = utils.decompress_pickle(solution)
                    seed_pop.update(pop)
                # Sort
                seed_pop = sorted(seed_pop)
                # Select the best and get only the RDKit molecule object
                MainConfig['seed_mol'] = [individual.mol for individual in seed_pop[:MainConfig['popsize']]]
            else:
                # Delete repeated SMILES
                MainConfig['seed_mol'] = set(MainConfig['seed_mol'])
                # COnvert to mol
                MainConfig['seed_mol'] = [Chem.MolFromSmiles(smi) for smi in MainConfig['seed_mol']]
                # Filter out invalid molecules
                MainConfig['seed_mol'] = list(filter(None, MainConfig['seed_mol']))
        else: # It will be assumed that it is a valid SMILES string
            MainConfig['seed_mol'] = Chem.MolFromSmiles(MainConfig['seed_mol'])

    # Convert if needed constraint_ref
    if 'constraint_ref' in MainConfig['costfunc_kwargs']:
        MainConfig['costfunc_kwargs']['constraint_ref'] = Chem.MolFromMolFile(MainConfig['costfunc_kwargs']['constraint_ref'])

    InitArgs = MainConfig.copy()

    # Modifying InitArgs
    _ = [InitArgs.pop(key, None) for key in ['type', 'njobs', 'pick']]
    InitArgs['costfunc'] = Cost

    # Getting call arguments
    CallArgs = dict()
    for key in ['njobs', 'pick']:
        try:
            CallArgs[key] = MainConfig[key]
        except KeyError:
            pass

    # Checking for follow jobs and sanity check on the arguments
    if FollowConfig:
        # Defining the possible mutable arguments with its default values depending on the type of run
        if MainConfig['type'].lower() == 'local':
            raise ValueError("Type = Local does not accept multiple call from the command line! Remove follow "\
                "jobs from the yaml file (only the main job is possible)")
        else:
            mutable_args = {
                'njobs': CallArgs['njobs'],
                'crem_db_path': InitArgs['crem_db_path'],
                'maxiter': InitArgs['maxiter'],
                'popsize': InitArgs['popsize'],
                'beta': InitArgs['beta'],
                'pc': InitArgs['pc'],
                'get_similar': InitArgs['get_similar'],
                # This one it will update with the default values of crem rather thant the previous one.
                'mutate_crem_kwargs': InitArgs['mutate_crem_kwargs'],
                'save_pop_every_gen': InitArgs['save_pop_every_gen'],
                'deffnm': InitArgs['deffnm'],
            }


        # Sanity check
        for job in FollowConfig:
            for arg in FollowConfig[job]:
                if arg not in mutable_args:
                    raise ValueError(f"The job: {job} has a non-valid argument \"{arg}\". "\
                        f"For now only the following are accepted: {list(mutable_args.keys())}")

    # Initialize the class
    ResultsClass = TypeOfRun(**InitArgs)
    # Call the class
    print(f"You are using moldrug: {__version__}.\n\nThe main job is being executed.\n\n"\
        f"Started at {datetime.datetime.now().strftime('%c')}")
    ResultsClass(**CallArgs)
    # Saving data
    if MainConfig['type'].lower() == 'local':
        ResultsClass.pickle("local_result", compress=True)
        utils.make_sdf(ResultsClass.pop, sdf_name = "local_pop")
    else:
        ResultsClass.pickle(f"{InitArgs['deffnm']}_result", compress=True)
        utils.make_sdf(ResultsClass.pop, sdf_name = f"{InitArgs['deffnm']}_pop")
    print('The main job finished!')
    # In case that follows jobs were defined
    if FollowConfig:
        for job in FollowConfig:
            print(f"The follow job {job} started.")

            # Updating arguments
            mutable_args.update(FollowConfig[job])
            InitArgs = mutable_args.copy()

            # Getting call arguments
            CallArgs = dict()
            for key in ['njobs', 'pick']:
                try:
                    CallArgs[key] = InitArgs[key]
                except KeyError:
                    pass

            # Changing the attributes values
            for arg in InitArgs:
                setattr(ResultsClass, arg, InitArgs[arg])

            # Call the class again
            ResultsClass(**CallArgs)
            # Saving data
            ResultsClass.pickle(f"{InitArgs['deffnm']}_result", compress=True)
            utils.make_sdf(ResultsClass.pop, sdf_name = f"{InitArgs['deffnm']}_pop")
            print(f'The job {job} finished!')

def __constraintconf_cmd():
    """
    Command line implementation for :meth:`moldrug.constraintconf.constraintconf`
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '--pdb',
        help = 'Protein pdb file',
        dest = 'pdb',
        type = str,
    )
    parser.add_argument(
        '--smi',
        help='Input SMILES file name',
        dest = 'smi',
        type = str,
    )
    parser.add_argument(
        '--fix',
        help = 'File with fixed piece of the molecule',
        dest = 'fix',
        type = str,
    )
    parser.add_argument(
        '--out',
        help = 'Output file name',
        dest = 'out',
        type = str,
    )
    parser.add_argument(
        '--max',
        help = 'Maximum number of conformers to generate, by default %(default)s',
        dest = 'max',
        default = 25,
        type = int,
    )
    parser.add_argument(
        '--rms',
        help = 'RMS cutoff, by default %(default)s',
        dest = 'rms',
        default = 0.01,
        type = float,
    )
    parser.add_argument(
        '--bump',
        help = 'Bump cutoff, by default %(default)s',
        dest = 'bump',
        default = 1.5,
        type = float,
    )
    args = parser.parse_args()
    constraintconf.constraintconf(args.pdb, args.smi, args.fix, args.out, args.max, args.rms, args.bump)


if __name__ == '__main__':
    pass