# model3_figure.py
#
# By Murillo F. Rodrigues, 25 August 2025
# Wall Lab, Oregon Health & Science University

from pathlib import Path
import tskit, os
import polars as pl
import numpy as np
import altair as alt

alt.data_transformers.enable("vegafusion")

# NOTE: the base repository path needs to be configured for your setup here!
repository_path = Path('/Users/rodrigmu/Documents/SimHumanity')

# set the current working directory to the SimHumanity repository
os.chdir(repository_path)

trees_path = repository_path / "simhumanity_trees_RO"

MU_TOTAL = 2.0e-8
WIN_LEN = 1_000_000

dfs = []

# Looping over chromosomes
for chrom in range(1, 23):
    chrom = str(chrom)
    print(f"Processing chromosome {chrom}...")

    # load the tree sequence
    print(f"   loading...")
    ts = tskit.load(trees_path / f"chromosome_{chrom}.trees")
    seq_len = ts.sequence_length

    # load the recombination map
    print(f"   loading recombination map...")
    rec_path = Path(f'stdpopsim extraction/extracted/chr{chrom}_recombination.txt')
    rec_df = pl.read_csv(rec_path, has_header=False, new_columns=["pos", "rate"])
    ends = np.append(rec_df["pos"].to_numpy()[1:], seq_len)
    rec_df = rec_df.with_columns(end = ends)

    # fetching windows with super low recombination rates to mask from figure
    low_rec_intervals = rec_df.filter(pl.col("rate") < 1e-12).select(["pos", "end"]).to_numpy()

    # masking the tree sequence at low recombination rate intervals, because of high variance in these regions
    print(f"   masking tree sequence...")
    ts = ts.delete_intervals(low_rec_intervals)


    # load the genomic elements
    print(f"   loading genomic elements...")
    ge_path = Path(f'stdpopsim extraction/extracted/chr{chrom}_genomic_elements.txt')
    ge_df = pl.read_csv(ge_path, has_header=False, new_columns=["type", "start", "end"])
    ge_df = ge_df.with_columns(pl.lit(chrom).alias("chrom"))

    # compute branch and site diversity within genomic elements
    print(f"   computing diversity...")
    breakpoints = np.append(ge_df["start"].to_numpy(), seq_len)
    site_div = ts.diversity(mode='site', windows=breakpoints)
    branch_div = ts.diversity(mode='branch', windows=breakpoints)
    ge_df = ge_df.with_columns(site_div = site_div, branch_div = branch_div, type = pl.when(pl.col("type") == 0).then(pl.lit("neutral")).otherwise(pl.lit("exon")).alias("type")).with_columns(ratio = pl.col("site_div")/pl.col("branch_div")) # create labels and compute the ratio of site/branch diversity
    dfs.append(ge_df)
    # saving the first chromosome for later use
    if chrom == "1":
        ts1 = ts
        ge_df1 = ge_df
        rec_df1 = rec_df

ge_df = pl.concat(dfs)
ge_df = ge_df.filter(pl.col("branch_div") > 0.01) # 0 branch diversity means low recombination
ge_df = ge_df.with_columns(log_ratio = (pl.col("ratio").log10()))

IN_TO_PX = 72
site_boxplot = alt.Chart(ge_df.with_columns(log_site_div = (pl.col("site_div")).log10())).mark_boxplot().encode(
    x = alt.X("type", title="", sort='descending'),
    y = alt.Y("log_site_div", title="Site diversity (log10-scaled)"),
    color = alt.Color("type", title="Genomic element type", sort="descending", legend=None),
).properties(
    width=11*IN_TO_PX*0.4,
    height=8.5*IN_TO_PX*0.4,
)
site_boxplot.save("site_boxplot.pdf")

branch_boxplot = alt.Chart(ge_df.filter(pl.col("branch_div") > 0.01).with_columns(log_branch_div = (pl.col("branch_div")).log10())).mark_boxplot().encode(
    x = alt.X("type", title="", sort='descending'),
    y = alt.Y("log_branch_div", title="Branch diversity (log10-scaled)", scale=alt.Scale(domain=[2, 6])),
    color = alt.Color("type", title="Genomic element type", sort="descending", legend=None),
).properties(
    width=11*IN_TO_PX*0.4,
    height=8.5*IN_TO_PX*0.4,
)
branch_boxplot.save("branch_boxplot.pdf")

# Plot of ratio site/branch div along chromosome 1
seq_len = ts1.sequence_length
breakpoints = np.arange(0,seq_len, step=WIN_LEN)
breakpoints = np.append(breakpoints, seq_len)

site_div = ts1.diversity(windows=breakpoints, mode="site")
branch_div = ts1.diversity(windows=breakpoints, mode="branch")

df = pl.DataFrame({"start":breakpoints[:-1], "end":breakpoints[1:], "site_div":site_div, "branch_div":branch_div})
df = df.with_columns(ratio = pl.col("site_div")/pl.col("branch_div"))

ratio_plot = alt.Chart(df.filter(pl.col("branch_div") > 0.01)).mark_line().encode(
    x = alt.X("start", title="Position"),
    y = alt.Y("ratio", title="Site/Branch Diversity Ratio", axis=alt.Axis(format=".1e"), scale=alt.Scale(domain=[1.8e-8, 2.2e-8])),
).properties(
    width=11*IN_TO_PX*0.8,
    height=8.5*IN_TO_PX*0.4,
)
ratio_plot.save("ratio_plot.pdf")

# Combined panel with site, branch and along chr1 plots
combined = alt.vconcat(site_boxplot.properties(title='A') | branch_boxplot.properties(title='B') , ratio_plot.properties(title='C') , center=True).configure_title(anchor='start')
combined.save("combined.png", ppi=400)
combined

