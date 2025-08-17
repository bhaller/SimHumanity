# model3_analysis.py
#
# By Ben Haller, 22 June 2025
# Messer Lab, Cornell University

# note that this requires pyslim 1.1.0 or later to run!

from pathlib import Path
import tskit, msprime, pyslim, warnings
from timeit import default_timer as timer


# NOTE: the base repository path needs to be configured for your setup here!
directory_path = Path('/Users/bhaller/Documents/Research/MesserLab/SLiM_project/Publication 2025 HumanPopGen/SimHumanity/simhumanity_trees')

# suppress time unit warnings from msprime; one generation equals one tick in this model, so it is fine
# ********** TO DO: instead of this, do timeUnits='generations' in the SLiM script so the tree sequence is marked as being in generations
warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

start = timer()

try:
    # iterate over the .trees files in the trees archive, in alphabetical order
    for file_path in sorted(directory_path.iterdir(), key=lambda x:x.name):
        if file_path.is_file() and file_path.suffix == '.trees':
            print(f"Processing {file_path.name}...")
            print(f"   loading...")
            ts = tskit.load(file_path)
            
            # remove vacant nodes from the sample, to avoid issues with "missing data" in stats computations
            ts = pyslim.remove_vacant(ts)
            
            # ********** TO DO: need to use the correct recombination rate map, as in the SLiM model!
            print(f"   recapitating...")
            ts = pyslim.recapitate(ts, ancestral_Ne=7310, recombination_rate=1e-8, random_seed=1)
            
            print(f"   computing diversity(mode='branch')...")
            d_branch = ts.diversity(mode='branch')
            
            # ********** TO DO: correct neutral mutation rates need to be used for each chromosome,
            # including the elevated mutation rate for the MT chromosome!
            print(f"   overlaying mutations...")
            next_id = pyslim.next_slim_mutation_id(ts)
            ts = msprime.sim_mutations(ts, rate=1e-7, random_seed=1,
                model=msprime.SLiMMutationModel(type=0, next_id=next_id), keep=True)
            
            print(f"   computing diversity(mode='site')...")
            d_site = ts.diversity(mode='site')
            
            # with correct mutation rates for the neutral mutation overlay we should see
            # that d_site = d_branch * mutation rate, roughly, although that is only the
            # exact expectation if all mutations are neutral
            print(f"       ", d_site, " ", d_branch)
except FileNotFoundError:
    print(f"Error: Directory '{directory_path}' not found.")
except Exception as e:
    print(f"An error occurred: {e}")

end = timer()
print("Elapsed time: ", end - start)
