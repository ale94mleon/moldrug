#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""For information of MolDrug:
    Docs: https://moldrug.readthedocs.io/en/latest/
    Source Code: https://github.com/ale94mleon/moldrug
"""
from moldrug._version import __version__, __version_tuple__
import yaml, argparse, inspect, importlib

def run():
    from moldrug import utils
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        help='The configuration yaml file',
        dest='yaml_file',
        type=str)
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('-f', '--fitness',
                        help="The path to the user-custom fitness module; inside of which the given custom cost function must be implemented. "\
                            "See the docs for how to do it properly. E.g. my/awesome/fitness_module.py."\
                            "By default will look in the moldrug.fitness module.",
                        dest='fitness',
                        default=None,
                        type=str)
    args = parser.parse_args()
    
    with open(args.yaml_file, 'r') as c:
        Config = yaml.safe_load(c)
    
    # Checking how many jobs are defined: Main and follows jobs
    MainConfig = Config.pop(list(Config.keys())[0])
    FollowConfig = Config

    if args.fitness:
        spec=importlib.util.spec_from_file_location('fitness', args.fitness)
        fitness = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(fitness)
        Cost = dict(inspect.getmembers(fitness))[MainConfig['costfunc']]
    else:
        from moldrug import fitness
        Cost = dict(inspect.getmembers(fitness))[MainConfig['costfunc']]
    
    if MainConfig['type'].lower() == 'ga':
        TypeOfRun = utils.GA
    elif MainConfig['type'].lower() == 'local':
        TypeOfRun = utils.Local
    else:
        raise RuntimeError(f"\"{MainConfig['type']}\" it is not a possible type. Select from: GA or Local")
    InitArgs = MainConfig.copy()

    # Modifying InitArgs
    [InitArgs.pop(key, None) for key in ['type', 'njobs']]
    InitArgs['costfunc'] = Cost

    # Checking for follow jobs and sanity check on the arguments
    if FollowConfig:
        # Defining the possible mutable arguments with its default values
        mutable_args = {
            'njobs': MainConfig['njobs'],
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
                    raise ValueError(f"The job: {job} has a non-valid argument \"{arg}\". For now only the following are accepted: {list(mutable_args.keys())}")   
    
    # Initialize the class
    ResultsClass = TypeOfRun(**InitArgs)
    # Call the class
    print(f'The main job is being executed.')
    ResultsClass(MainConfig['njobs'])
    # Save final data
    ResultsClass.pickle(f"{InitArgs['deffnm']}_result", compress=True)
    # Saving final sdf file always
    utils.make_sdf(ResultsClass.pop, sdf = f"{InitArgs['deffnm']}_pop.sdf")
    print(f'The main job finished!.')
    # In case that follows jobs were defined
    if FollowConfig:
        for job in FollowConfig:
            print(f"The follow job {job} started.")
            
            # Updating arguments
            mutable_args.update(FollowConfig[job])
            InitArgs = mutable_args.copy()
            njobs = InitArgs.pop('njobs')

            # Changing the attributes values
            for arg in InitArgs:
                setattr(ResultsClass, arg, InitArgs[arg])

            # Call the class again
            ResultsClass(njobs)
            # Save final data
            ResultsClass.pickle(f"{InitArgs['deffnm']}_result", compress=True)
            # Saving final sdf file always
            utils.make_sdf(ResultsClass.pop, sdf = f"{InitArgs['deffnm']}_pop.sdf")
            print(f'The job {job} finished!.')