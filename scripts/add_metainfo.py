import pandas as pd

table = pd.read_csv(snakemake.input[0], sep='\t')

genetypes = pd.read_csv(snakemake.params["genetype"], sep='\t')

table = table.merge(genetypes, on="gene", how="left")

table.to_csv(snakemake.output[0], sep='\t', index=False)
