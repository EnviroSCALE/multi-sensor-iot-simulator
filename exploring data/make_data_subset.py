
import pandas as pd
import numpy as np


#0	  1		2		3		4	5		6		7			8	9
#Date,Time,device_id,PM2.5,PM10,PM1,Temperature,Humidity, lat, lon
#2017-01-01,08:00:00,74DA388FF60A,31,33,22,22.75,78,25.072,121.657


df = pd.read_csv('201701_Taiwan.csv',  delimiter=',', header=0)
print "reading done"
####Count_Row=df.shape[0] #gives number of row count
#####Count_Col=df.shape[1] #gives number of col count
#####print Count_Row
#filter rows
df=df.iloc[0:1000000]
df.to_csv("shorted_humidity_1M.csv", columns=['Humidity'], header=True)
