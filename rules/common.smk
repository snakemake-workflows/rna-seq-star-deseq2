def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def exist_strandedness(units):
    return "strandedness" in units.columns
    
def strandedness(units):
    strandedness_list = []
    if exist_strandedness(units):
        for unit in units.itertuples():
            strand_val = units.loc[(unit.sample, unit.unit), "strandedness"]
            if pd.isnull(strand_val) or strand_val == "0":
                strandedness_list.append(1) #non stranded protocol
            elif strand_val == "yes":
                strandedness_list.append(2) #3rd column
            elif strand_val == "reverse":
                strandedness_list.append(3) #4th column, usually for Illumina truseq
            else:
                raise ValueError(("'strandedness' column should be empty or have the " 
                "value 0, 'yes' or 'reverse',instead has the value {}").format(repr(strand_val)))
    else:
        strandedness_list.append(1)#non stranded for cases where there isn't a "strandedness" column
        strandedness_list *= units.shape[0]
    
    return strandedness_list
