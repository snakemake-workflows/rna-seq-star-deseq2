import pandas as pd

def strandedness(strand_vals):
    strandedness_list = []
    for strand_val in snakemake.params.strand:
        if pd.isnull(strand_val) or strand_val == "none":
            strandedness_list.append(1) #non stranded protocol
        elif strand_val == "yes":
            strandedness_list.append(2) #3rd column
        elif strand_val == "reverse":
            strandedness_list.append(3) #4th column, usually for Illumina truseq
        else:
            raise ValueError(("'strandedness' column should be empty or have the " 
                "value 'none', 'yes' or 'reverse', instead has the value {}").format(repr(strand_val)))
    
    return strandedness_list

stranded_list = strandedness(snakemake.params.strand)
counts = [pd.read_table(f, index_col=0, usecols=[0, coln], header=None, skiprows=4) \
for f, coln in zip(snakemake.input, stranded_list)]

for t, sample in zip(counts, snakemake.params.samples):
    t.columns = [sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
# collapse technical replicates
matrix = matrix.groupby(matrix.columns, axis=1).sum()
matrix.to_csv(snakemake.output[0], sep="\t")
