import pandas as pd

counts = [pd.read_table(f, index_col=0, usecols=[0, 1], header=None, skiprows=4)
          for f in snakemake.input]

for t, (sample, unit) in zip(counts, snakemake.params.units.index):
    t.columns = [sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
# collapse technical replicates
matrix = matrix.groupby(matrix.columns, axis=1).sum()
matrix.to_csv(snakemake.output[0], sep="\t")
