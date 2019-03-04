import pandas as pd

counts = [pd.read_table(f, index_col=0, usecols=[0, 1], header=None, skiprows=4)
          for f in snakemake.input]

for t, unit in zip(counts, units):
    t.columns = [unit.sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
# collapse technical replicates
matrix = matrix.groupby(matrix.columns, axis=1).sum()
print(matrix)
matrix.to_csv(snakemake.output[0], sep="\t")
