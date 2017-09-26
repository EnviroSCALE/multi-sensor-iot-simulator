
import pandas as pd

import numpy as np
import sys

#0	  1		2		3		4	5		6		7			8	9
#Date,Time,device_id,PM2.5,PM10,PM1,Temperature,Humidity, lat, lon
#2017-01-01,08:00:00,74DA388FF60A,31,33,22,22.75,78,25.072,121.657
#input_data = "hist_data.csv"
#output_png = "histogram.png"	
#t = np.genfromtxt('/home/tasnim/Desktop/201701_Taiwan.csv', delimiter=',')
#t2 = np.genfromtxt('cv.csv', delimiter=',')
#df = pd.read_csv('/home/tasnim/Desktop/201701_Taiwan.csv',  delimiter=',', header=0)
#t = df['Humidity'].values

'''
df = pd.read_csv('/home/tasnim/Desktop/201701_Taiwan.csv',  delimiter=',', header=0)
print "read done"
####Count_Row=df.shape[0] #gives number of row count
#####Count_Col=df.shape[1] #gives number of col count
#####print Count_Row
#filter rows
df=df.iloc[0:1000000]
df.to_csv("shorted_dust1M.csv", columns=['PM2.5'], header=True)
'''


df = pd.read_csv('shorted_dust1M.csv',  delimiter=',', header=0)
#df = pd.read_csv('aqdata.csv',  delimiter=',', header=0)
temp = df['PM2.5'].values
print temp
#temp= t[:, 6]
#temp = temp[0:25]
#print temp

#import matplotlib.pyplot as plt
import pylab as plt
import numpy as np
import matplotlib.mlab as mlab

#Temp = t[:, 6]
#D = D[0:25]
#TotError = TotError[0:25]

#hist = np.histogram(temp, bins=np.arange(np.min(temp), np.max(temp), 1), density=True)
#n, bins, patches = hist
#print hist
#n, bins, patches = plt.hist(temp, 50, normed=1, facecolor='blue', alpha=0.35)
weights = np.ones_like(temp)/float(len(temp))
n, bins, patches = plt.hist(temp, 500, edgecolor='black', weights=weights)

mu = np.mean(temp)
#sigma = np.power(np.var(temp),0.09) 
sigma = np.power(np.var(temp),0.03) 

# add a 'best fit' line
y = mlab.normpdf( bins, mu, sigma)
print n, bins
#y= 0.8*np.exp(-bins/np.mean(bins)) 
l = plt.plot(bins, y, 'r--', linewidth=5)


plt.xlabel('Sensor Reading')
plt.ylabel('Probability')
plt.title('Dust Histogram')
#plt.axis([np.min(temp), np.max(temp), 0, 1])
#plt.grid(True)
plt.savefig('dust_histogram_taiwan.pdf', format='pdf')
plt.show()
