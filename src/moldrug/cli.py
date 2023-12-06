#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For information of MolDrug:
    Docs: https://moldrug.readthedocs.io/en/latest/
    Source Code: https://github.com/ale94mleon/moldrug
"""
from moldrug import utils, constraintconf, __version__
import yaml
import argparse
import inspect
import os
import sys
import datetime
from rdkit import Chem
from typing import Union


class CommandLineHelper:
    def __init__(self, parser) -> None:
        self.args = parser.parse_args()
        self.yaml_file = self.args.yaml_file
        self.fitness = self.args.fitness
        self.continuation = self.args.continuation
        self.verbose = self.args.verbose
        self._set_attributes()

        if self.verbose:
            os.environ['MOLDRUG_VERBOSE'] = 'True'

    def _set_attributes(self):
        # Get and set configuration
        self._set_config()
        # It will generate costfunc attribute
        self._set_costfunc()
        # it will generate the attribute TypeOfRun
        self._set_TypeOfRun()
        # it will generate the attributes: MainConfig, FollowConfig, InitArgs, CallArgs and MutableArgs
        self._translate_config()
        # Here FollowConfig is updated and the corresponded initialization occurs if self.continuation.
        # pbz2 and new_maxiter attributes are generated
        self._set_init_MolDrugClass()

    def _set_config(self):
        with open(self.yaml_file, 'r') as c:
            self.configuration = yaml.safe_load(c)

    def _split_config(self):
        config = self.configuration.copy()
        MainConfig = config.pop(list(config.keys())[0])
        FollowConfig = config
        return MainConfig, FollowConfig

    def _set_costfunc(self):
        if self.fitness:
            # If the fitness module provided is not in the current directory or if its name is not fitness
            with open(self.fitness, 'r') as source:
                with open('CustomMolDrugFitness.py', 'w') as destination:
                    destination.write(source.read())
            sys.path.append('.')
            import CustomMolDrugFitness
            costfunc = dict(inspect.getmembers(CustomMolDrugFitness))[self._split_config()[0]['costfunc']]
        else:
            from moldrug import fitness
            costfunc = dict(inspect.getmembers(fitness))[self._split_config()[0]['costfunc']]
        self.costfunc = costfunc

    def _set_TypeOfRun(self):
        self._TypeOfRun_str = self._split_config()[0]['type'].lower()
        if self._TypeOfRun_str == 'ga':
            self.TypeOfRun = utils.GA
        elif self._TypeOfRun_str == 'local':
            self.TypeOfRun = utils.Local
        else:
            raise NotImplementedError(f"\"{self._split_config()[0]['type']}\" it is not a possible type. "
                                      "Select from: GA or Local")

    def _translate_config(self):
        MainConfig, FollowConfig = self._split_config()

        # Convert the SMILES (or path to compress_pickle) to RDKit mol (or list of RDkit mol)
        if self._TypeOfRun_str == 'local':
            MainConfig['seed_mol'] = Chem.MolFromSmiles(MainConfig['seed_mol'])
        else:
            # TODO make clear errors in case the path does not exist or it was not possible to create a molecule
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
            else:  # It will be assumed that it is a valid SMILES string
                MainConfig['seed_mol'] = Chem.MolFromSmiles(MainConfig['seed_mol'])

        # Convert if needed constraint_ref
        if 'constraint_ref' in MainConfig['costfunc_kwargs']:
            MainConfig['costfunc_kwargs']['constraint_ref'] = Chem.MolFromMolFile(
                MainConfig['costfunc_kwargs']['constraint_ref'])

        InitArgs = MainConfig.copy()

        # Modifying InitArgs
        _ = [InitArgs.pop(key, None) for key in ['type', 'njobs', 'pick']]
        InitArgs['costfunc'] = self.costfunc

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
                raise ValueError("Type = Local does not accept multiple call from the command line! Remove follow "
                                 "jobs from the yaml file (only the main job is possible)")
            else:
                # Add default value in case it is not provided for keyword arguments
                list_of_keywords = [
                    'beta', 'pc', 'get_similar', 'mutate_crem_kwargs',
                    'save_pop_every_gen', 'checkpoint', 'deffnm',
                ]
                for param in inspect.signature(self.TypeOfRun).parameters.values():
                    if (param.kind == param.POSITIONAL_OR_KEYWORD and
                            param.default is not param.empty and
                            param.name in list_of_keywords and
                            param.name not in InitArgs):
                        InitArgs[param.name] = param.default

                MutableArgs = {
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
                    'checkpoint': InitArgs['checkpoint'],
                    'deffnm': InitArgs['deffnm'],
                }

            # Sanity check
            for job in FollowConfig:
                for arg in FollowConfig[job]:
                    if arg not in MutableArgs:
                        raise ValueError(f"The job: {job} has a non-valid argument \"{arg}\". "
                                         f"For now only the following are accepted: {list(MutableArgs.keys())}")
        else:
            MutableArgs = None

        self.MainConfig = MainConfig
        self.FollowConfig = FollowConfig
        self.InitArgs = InitArgs
        self.CallArgs = CallArgs
        self.MutableArgs = MutableArgs

    def _get_continuation_point(self):
        """
        This gave me the job and how many generation are needed to complete it.
        The further jobs are suppose that must run.
        """
        if self.continuation:
            if self._TypeOfRun_str != 'ga':
                raise RuntimeError('Continuation is only valid for GA runs.')
            # Check what was already done
            total_iter = 0
            pbz2 = None
            for job in self.configuration:
                # El problema es que no estan los archivos entonces hay que modificar total_iter
                if os.path.isfile(f"{self.configuration[job]['deffnm']}_result.pbz2"):
                    # must be defined maxiter in the configuration file
                    total_iter += self.configuration[job]['maxiter']
                    # Delete (update) the jobs in self.FollowConfig if they were already done
                    if job in self.FollowConfig:
                        del self.FollowConfig[job]
                    # Stay with the last one
                    pbz2 = f"{self.configuration[job]['deffnm']}_result.pbz2"

            # If there is a continuation file, use this
            if os.path.isfile("cpt.pbz2"):
                pbz2 = 'cpt.pbz2'
                iter_done = utils.decompress_pickle(pbz2).NumGens
                total_iter = 0
                for job in self.configuration:
                    total_iter += self.configuration[job]['maxiter']
                    if total_iter >= iter_done:
                        del self.FollowConfig[job]
                        break
            elif pbz2:
                iter_done = utils.decompress_pickle(pbz2).NumGens
            else:
                iter_done = 0
            new_maxiter = total_iter - iter_done
        else:
            # In this case we must start from scratch. There are not .pbz2 files in the directory
            pbz2, new_maxiter = None, 0

        # Set the attributes
        self.pbz2 = pbz2
        self.new_maxiter = new_maxiter

    def _set_init_MolDrugClass(self):
        # Here is where the continuation code is added

        # Get if if needed to continue and make the corresponded updates on self.FollowConfig
        self._get_continuation_point()

        if self.pbz2:
            self.MolDrugClass = utils.decompress_pickle(self.pbz2)
            self.MolDrugClass.maxiter = self.new_maxiter
        else:
            # Initialize the class from scratch
            self.MolDrugClass = self.TypeOfRun(**self.InitArgs)

    def run_MolDrugClass(self):
        self.MolDrugClass(**self.CallArgs)

    def save_data(self):
        # Saving data
        if self._TypeOfRun_str == 'local':
            self.MolDrugClass.pickle("local_result", compress=True)
            utils.make_sdf(self.MolDrugClass.pop, sdf_name="local_pop")
        else:
            self.MolDrugClass.pickle(f"{self.MolDrugClass.deffnm}_result", compress=True)
            utils.make_sdf(self.MolDrugClass.pop, sdf_name=f"{self.MolDrugClass.deffnm}_pop")

    def __repr__(self) -> str:
        string = self.args.__repr__().replace('Namespace', self.__class__.__name__)
        if self.continuation:
            string += f"\nContinuationPoint(pbz2={self.pbz2}, do_iter={self.new_maxiter})"
        return string


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
    parser.add_argument("-f", "--fitness",
                        help="The path to the user-custom fitness module; inside of which the given custom "
                        "cost function must be implemented. "
                        "See the docs for how to do it properly. E.g. my/awesome/fitness_module.py. "
                        "By default will look in the moldrug.fitness module.",
                        dest="fitness",
                        nargs=argparse.OPTIONAL,
                        default=None,
                        type=str)
    parser.add_argument("-c", "--continue",
                        help="To continue the simulation. The MolDrug command must be the same "
                        "and all the output MolDrug files must be located "
                        "in the working directory. This option is only compatible "
                        "with moldrug.utils.GA; otherwise, a RuntimeError will be raised.",
                        action="store_true",
                        dest="continuation")
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f"moldrug: {__version__}")
    parser.add_argument(
        '-V', '--verbose',
        nargs="?",
        dest='verbose',
        const=True,
        default=False,
        type=bool
    )
    UserArgs = CommandLineHelper(parser)

    print(
        f"Started at {datetime.datetime.now().strftime('%c')}\n"
        f"You are using moldrug: {__version__}.\n\n"
        f"{UserArgs}\n\n")

    # Call the class
    UserArgs.run_MolDrugClass()
    # Saving data
    UserArgs.save_data()
    # print('The main job finished!')

    # In case that follows jobs were defined
    if UserArgs.FollowConfig:
        MutableArgs = UserArgs.MutableArgs.copy()
        for job in UserArgs.FollowConfig:
            print(f"The follow job {job} started.")

            # Updating arguments
            MutableArgs.update(UserArgs.FollowConfig[job])
            InitArgs = MutableArgs.copy()

            # Changing the attributes values
            for arg in InitArgs:
                setattr(UserArgs.MolDrugClass, arg, InitArgs[arg])

            # Call the class again
            UserArgs.run_MolDrugClass()
            # Saving data
            UserArgs.save_data()
            print(f'The job {job} finished!')

    # Clean checkpoint on normal end
    if os.path.isfile('cpt.pbz2'):
        os.remove('cpt.pbz2')


def __constraintconf_cmd():
    """
    Command line implementation for :meth:`moldrug.constraintconf.constraintconf`
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--pdb",
        help="Protein pdb file",
        dest="pdb",
        type=str,
    )
    parser.add_argument(
        "--smi",
        help="Input SMILES file name",
        dest="smi",
        type=str,
    )
    parser.add_argument(
        "--fix",
        help="File with fixed piece of the molecule",
        dest="fix",
        type=str,
    )
    parser.add_argument(
        "--out",
        help="Output file name",
        dest="out",
        type=str,
    )
    parser.add_argument(
        "--max",
        help="Maximum number of conformers to generate, by default %(default)s",
        dest="max",
        default=25,
        type=int,
    )
    parser.add_argument(
        "--rms",
        help="RMS cutoff, by default %(default)s",
        dest="rms",
        default=0.01,
        type=float,
    )
    parser.add_argument(
        "--bump",
        help="Bump cutoff, by default %(default)s",
        dest="bump",
        default=1.5,
        type=float,
    )
    parser.add_argument(
        "--seed",
        help="Provide a seed for the random number generator so that "
        "the same coordinates can be obtained for a molecule on multiple runs. "
        "If None, the RNG will not be seeded, by default None %(default)s",
        dest="seed",
        default=None,
        type=Union[int, None],
    )
    args = parser.parse_args()
    constraintconf.constraintconf(
        pdb=args.pdb,
        smi=args.smi,
        fix=args.fix,
        out=args.out,
        max_conf=args.max,
        rms=args.rms,
        bump=args.bump,
        randomseed=args.seed)


if __name__ == '__main__':
    pass
