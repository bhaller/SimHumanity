# model3_analysis.py
#
# By Ben Haller and Murillo F. Rodrigues, 24 August 2025
# Messer Lab, Cornell University

# note that this requires pyslim 1.1.0 or later to run!

from pathlib import Path
import tskit, msprime, pyslim, warnings, os
from timeit import default_timer as timer
import polars as pl
import numpy as np

# NOTE: the base repository path needs to be configured for your setup here!
repository_path = Path('/Users/bhaller/Documents/Research/MesserLab/SLiM_project/Publication 2025 HumanPopGen/SimHumanity')

# set the current working directory to the SimHumanity repository
os.chdir(repository_path)

out_path = repository_path / "simhumanity_trees_RO"
out_path.mkdir(parents=False, exist_ok=False)

MU_TOTAL = 2.0e-8
MU_BENEF = 1e-12
MU_DELET = 1.2e-8

def get_recomb_map(chrom, seq_len):
    """
    Returns a recombination rate map (msprime.RateMap) for a given chromosome.
    Note the relative path to stdpopsim extraction is hardcoded in the function.
    """
    if chrom == 'MT':
        positions = np.array([0, seq_len])
        rates = np.array([0.0])
    else:
        rec_path = Path(f'stdpopsim extraction/extracted/chr{chrom}_recombination.txt')
        rec_df = pl.read_csv(rec_path, has_header=False, new_columns=["pos", "rate"])

        positions = np.append(rec_df['pos'].to_numpy(), seq_len) # need to include the end of the sequence
        rates = rec_df['rate'].to_numpy()
    return msprime.RateMap(position=positions, rate=rates)

def get_mut_map(chrom, seq_len, mt_multiplier = 20):
    """
    Returns a mutation rate map (msprime.RateMap) for a given chromosome.
    For the MT chromosome, the mutation rate is constant and equal to MU_TOTAL multiplied by mt_multiplier.
    Note the relative paths to stdpopsim extraction and the mutation rates are hardcoded in the function.
    """
    if chrom == 'MT':
        breakpoints = np.array([0, seq_len])
        rates = np.array([MU_TOTAL*mt_multiplier])
    else:
        ge_path = Path(f'stdpopsim extraction/extracted/chr{chrom}_genomic_elements.txt')
        ge_df = pl.read_csv(ge_path, has_header=False, new_columns=["type", "start", "end"])
        ge_df = ge_df.with_columns(pl.when(pl.col("type") == 0).then(pl.lit(MU_TOTAL)).otherwise(pl.lit(MU_TOTAL-MU_BENEF-MU_DELET)).alias("neutral_rate")) # computing the neutral rate for each genomic element
        breakpoints = np.append(ge_df["start"].to_numpy(), seq_len)
        rates = ge_df["neutral_rate"].to_numpy()
    return msprime.RateMap(position=breakpoints, rate=rates)


start = timer()

try:
    # iterate over the .trees files in the trees archive, in alphabetical order
    for file_path in sorted((repository_path / "simhumanity_trees").iterdir(), key=lambda x:x.name):
        if file_path.is_file() and file_path.suffix == '.trees':
            print(f"Processing {file_path.name}...")
            print(f"   loading...")

            # extract the chromosome number from the file path
            chrom = str(file_path).split('.')[0].split('_')[-1]
            ts = tskit.load(file_path)

            # remove vacant nodes from the sample, to avoid issues with "missing data" in stats computations
            ts = pyslim.remove_vacant(ts)

            print(f"Basic tree sequence information:")
            print(ts)


            # recapitate; note this issues a warning about large provenance, because it includes the rate map in the provenance
            print(f"   fetching the recombination rate map...")
            rec_map = get_recomb_map(chrom, ts.sequence_length)

            print(f"   recapitating...")
            ts = pyslim.recapitate(ts, ancestral_Ne=7310, recombination_rate=rec_map, random_seed=1)


            # compute where the recapitation period starts (in time ago)
            recap_start_time = ts.metadata["SLiM"]["tick"]-ts.metadata["SLiM"]["cycle"]
            print(f"   recapitation period starts at {recap_start_time} time ago")

            num_slim_muts = ts.num_mutations
            print(f"   {num_slim_muts} SLiM mutations")


            # overlay neutral mutations during the recapitation period
            print(f"   overlaying mutations over the recapitated portion...")
            next_id = pyslim.next_slim_mutation_id(ts)

            ts = msprime.sim_mutations(ts, rate=MU_TOTAL, random_seed=1, model=msprime.SLiMMutationModel(type=0, next_id=next_id), keep=True, start_time=recap_start_time) # start_time is in time ago, so we need to add mutations with a constant rate across the chromosome starting at the time the recapitated period starts (and older)

            num_postrecap_muts = ts.num_mutations
            print(f"   {num_postrecap_muts-num_slim_muts} mutations overlaid over the recapitated portion")


            # overlay neutral mutations during the SLiM period
            print(f"   fetching the mutation rate map for the SLiM period...")
            mut_map = get_mut_map(chrom, ts.sequence_length)

            print(f"   overlaying mutations over the SLiM period...")
            next_id = pyslim.next_slim_mutation_id(ts)

            ts = msprime.sim_mutations(ts, rate=mut_map, random_seed=1, model=msprime.SLiMMutationModel(type=0, next_id=next_id), keep=True, end_time=recap_start_time)


            num_total_muts = ts.num_mutations
            print(f"   {num_total_muts-num_slim_muts} mutations overlaid over the SLiM portion for a total of {num_total_muts} mutations")


            # write out the modified tree sequence to simhumanity_trees_RO
            print(f"   writing modified tree sequence...")
            ts.dump(str(out_path / file_path.name))
except FileNotFoundError:
    print(f"Error: Directory '{directory_path}' not found.")
except Exception as e:
    print(f"An error occurred: {e}")

end = timer()
print("Elapsed time: ", end - start)
