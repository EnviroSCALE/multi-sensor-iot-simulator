import pandas as pd
import numpy as np


def filterco2():
	df = pd.read_csv('data_co22.csv', delimiter=',', header=0)
	print "Read complete"
	Count_Row=df.shape[0] #gives number of row count
	Count_Col=df.shape[1] #gives number of col count
	# filter rows
	#df=df.iloc[0:1000000]

	old_names = ['', 'PM2.5'] 
	new_names = ['', '1']
	df.rename(columns=dict(zip(old_names, new_names)), inplace=True)

	df.to_csv("data_co2.csv",  columns=['1'], index=False, header=True)



df = pd.read_csv('data1M.csv', delimiter=',', header=0)
print "Read complete"
Count_Row=df.shape[0] #gives number of row count
Count_Col=df.shape[1] #gives number of col count
# filter rows
#df=df.iloc[0:1000000]

old_names = ['', 'Humidity'] 
new_names = ['', '0.25']
df.rename(columns=dict(zip(old_names, new_names)), inplace=True)

df.to_csv("data_humidity_taiwan.csv",  columns=['0.25'], index=False, header=True)
