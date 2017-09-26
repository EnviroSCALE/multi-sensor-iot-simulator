#import matplotlib.pyplot as plt
import pylab as plt
import numpy as np
# Equal weights approach 1
# Accumulated CV approach  2
# Separated CV approach  3
# Separated CV approach  4
# lnx approach  5
# sqrt approach  6
# no approach  7

class SimRun:
	def __init__(self, name, input):
		self.input = input
		self.name = name

run = []
run.append(SimRun('6-UP 5', '6-uperiod5.csv'))
run.append(SimRun('6-UP 50', '6-uperiod50.csv'))
run.append(SimRun('7-UP 50', '7-uperiod50.csv'))
run.append(SimRun('1-UP 50', '1-uperiod50.csv'))
run.append(SimRun('5-UP 50', '5-uperiod50.csv'))
run.append(SimRun('2-UP 50', '2-uperiod50.csv'))
print(len(run))


t = [None] * len(run)
TotError = [None] * len(run)
Util = [None] * len(run)
D = [None] * len(run)
AvgDelay = [None] * len(run)

#self.tot, self.prominent_error, self.avg_delay, self.utilization, choice, time (s), data (bytes)
#0          1                       2               3               4           5       6

for i in range(len(run)):
	t[i] = np.genfromtxt(run[i].input, delimiter=',')
	TotError[i] = t[i][:, 0]
	Util[i] = t[i][:, 3]
	D[i] = t[i][:, 6]/1000000.0 #MB
	AvgDelay[i] = t[i][:, 2]

# preparing marker positions - drop every alternate marker
xx = (np.arange(0, len(D[0]), 2))
xx =(xx.astype(int))
xx=list(xx.tolist())
MARKERS = ['-bo', '-gD', '-rs', '-m*',  '-cp', '-yH', '-yH', '-rp', '-bx']

for i in range(len(run)):
	plt.plot(D[i], AvgDelay[i],  MARKERS[i], markevery=xx,  label=run[i].name)
plt.legend(loc='upper right')
plt.xlabel('Data (MB)')
plt.ylabel('delay')
plt.title('D vs delay, 10 hour run, 7:3 weight')
plt.grid(True)
plt.savefig("d vs delay.pdf", format='pdf')
plt.show()

for i in range(len(run)):
	plt.plot(D[i], TotError[i],  MARKERS[i], markevery=xx,  label=run[i].name)
plt.legend(loc='upper right')
plt.xlabel('Data (MB)')
plt.ylabel('E')
plt.title('D vs E, 10 hour run, 7:3 weight')
plt.grid(True)
plt.savefig("d vs e.pdf", format='pdf')
plt.show()

for i in range(len(run)):
	plt.plot(D[i], Util[i],  MARKERS[i], markevery=xx,  label=run[i].name)
plt.legend(loc='upper right')
plt.xlabel('Data (MB)')
plt.ylabel('Util')
plt.title('D vs Util, 10 hour run, 7:3 weight')
plt.grid(True)
plt.savefig("d vs util_cv.pdf", format='pdf')
plt.show()



#~ plt.plot(D, Util,  '-bo', markevery=xx,  label=run1.name)
#~ plt.plot(D, Util2,  '-gD', markevery= xx,  label=run2.name)
#~ plt.plot(D, Util3,  '-rs', markevery= xx,  label=run3.name)
#~ plt.legend(loc='upper right')
#~ plt.xlabel('Data (MB)')
#~ plt.ylabel('Util')
#~ plt.title('D vs Util, 10 hour run, 7:3 weight')
#~ plt.grid(True)
#~ plt.savefig("d vs util_cv.pdf", format='pdf')
#~ plt.show()
