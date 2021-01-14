def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])
