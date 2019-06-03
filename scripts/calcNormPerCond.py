import pandas as pd
from functools import reduce
import os



normCounts = pd.read_csv(snakemake.input[0],sep=',',index_col=0)
samples = pd.read_csv(snakemake.params['samples'],sep='\t')

result = pd.DataFrame(index=normCounts.index)


conditionsCounts = {}

for index, row in samples.iterrows():

	sampleName = row[0]
	
	conditions = '_'.join(row[1:]) #Simply use a concatenated string as an identifier for the condition
	
	#Keep track of the samples per condition to calculate mean later
	if not conditions in conditionsCounts:
		conditionsCounts[conditions] = 1
	else:
		conditionsCounts[conditions] += 1
		
	#Create column if it doesn't exist yet, else sum
	
	if not conditions in result:
		result[conditions] = normCounts[sampleName]
	else:
		result[conditions] += normCounts[sampleName]
	
#Means

for column in result:
	if column in conditionsCounts:
		result[column] = result[column].apply(lambda x: x/conditionsCounts[column])
	else:
		print(column,' has no associated samples')
		
result.to_csv(snakemake.output[0])

