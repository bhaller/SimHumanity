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

def violin_and_boxplot(data, y_col, x_col, title, IN_TO_PX=72):
    violins = alt.Chart().transform_density(
        y_col,
        as_=[y_col, 'density'],
        groupby=['type']
    ).mark_area(orient='horizontal').encode(
        y=alt.Y(y_col, title=title),
        color=alt.Color(x_col, legend=None),
        x=alt.X(
            'density:Q',
            stack='center',
            impute=None,
            title=None,
            axis=alt.Axis(labels=False, values=[0],grid=False, ticks=True),
            scale=alt.Scale(nice=False,zero=False),
        ),
    )

    boxplot = alt.Chart().mark_boxplot(size=5, extent=0, outliers=False).encode(
            y=alt.Y(y_col, title=title),
            color=alt.value('black')
        )

    violin_boxplot = alt.layer(
        violins,
        boxplot
    ).properties(
        width=11*IN_TO_PX*0.2,
        height=8.5*IN_TO_PX*0.4,

    ).facet(
        data=data,
        column=alt.Column(
            x_col,
            header=alt.Header(
                titleOrient='bottom',
                labelOrient='bottom',
                labelPadding=0,
            ),
            title=None
        )
    ).resolve_scale(x=alt.ResolveMode("independent"))
    return violin_boxplot

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

# Plot violins for site and branch diversity
data = ge_df.with_columns(
    log_site_div = pl.max_horizontal(pl.lit(1e-8), pl.col("site_div")).log10()
)
site_violin = violin_and_boxplot(data, "log_site_div", "type", "Site diversity (log10-scaled)")
site_violin.save("site_plot.pdf")

data=ge_df.with_columns(log_branch_div = (pl.col("branch_div")).log10())
branch_violin = violin_and_boxplot(data, "log_branch_div", "type", "Branch diversity (log10-scaled)")
branch_violin.save("branch_plot.pdf")
branch_violin

# Plot of ratio site/branch div along chromosome 1
seq_len = ts1.sequence_length
breakpoints = np.arange(0,seq_len, step=WIN_LEN)
breakpoints = np.append(breakpoints, seq_len)

site_div = ts1.diversity(windows=breakpoints, mode="site")
branch_div = ts1.diversity(windows=breakpoints, mode="branch")

df = pl.DataFrame({"start":breakpoints[:-1], "end":breakpoints[1:], "site_div":site_div, "branch_div":branch_div})
df = df.with_columns(ratio = pl.col("site_div")/pl.col("branch_div"))

ratio_plot = alt.Chart(df.filter(pl.col("branch_div") > 0.01)).mark_line(color="black").encode(
    x = alt.X("start", title="Position"),
    y = alt.Y("ratio", title="Site/Branch Diversity Ratio", axis=alt.Axis(format=".1e"), scale=alt.Scale(domain=[1.8e-8, 2.2e-8])),
).properties(
    width=11*IN_TO_PX*0.8,
    height=8.5*IN_TO_PX*0.4,
)
ratio_plot.save("ratio_plot.pdf")

# Combined panel with site, branch and along chr1 plots
combined = alt.vconcat(site_violin.properties(title='A') | branch_violin.properties(title='B') , ratio_plot.properties(title='C') , center=True).configure_title(anchor='start').configure_facet(
    spacing=0
).configure_view(
    stroke=None
)
combined.save("figure3.png", ppi=400)
combined.save("figure3.pdf")
