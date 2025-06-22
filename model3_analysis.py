// model3_analysis.py.slim
//
// By Ben Haller, 22 June 2025
// Messer Lab, Cornell University

from pathlib import Path
import tskit, msprime, pyslim


// NOTE: the base repository path needs to be configured for your setup here!
directory_path = Path('/Users/bhaller/Documents/Research/MesserLab/SLiM_project/Publication 2025 HumanPopGen/SimHumanity/simhumanity_trees')

try:
    # iterate over the .trees files in the trees archive, in alphabetical order
    for file_path in sorted(directory_path.iterdir(), key=lambda x:x.name):
        if file_path.is_file() and file_path.suffix == '.trees':
            print(f"Processing {file_path.name}...")
            print(f"   loading...")
            ts = tskit.load(file_path)
            
            # this errors, which I think is due to a bug in pyslim; stay tuned
            print(f"   recapitating...")
            ts = pyslim.recapitate(ts, ancestral_Ne=7310, recombination_rate=1e-8, random_seed=1)
            
            print(f"   computing diversity(mode='branch')...")
            d_branch = ts.diversity(mode='branch')
            
            print(f"   overlaying mutations...")
            ts = msprime.sim_mutations(ts, rate=1e-7, random_seed=1, keep=True)
            
            print(f"   computing diversity(mode='site')...")
            d_site = ts.diversity(mode='site')
            
            print(f"       ", d_site, " ", d_branch)
except FileNotFoundError:
    print(f"Error: Directory '{directory_path}' not found.")
except Exception as e:
    print(f"An error occurred: {e}")
