import pandas as pd

matrix = pd.concat([pd.read_table(f, index_col=0)[1] for f in snakemake.input], 
                   axis=1, names=snakemake.params.samples)
matrix.to_csv(snakemake.output[0], sep="\t")
