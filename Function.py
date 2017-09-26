import math
import numpy as np;

class TraceDataGenerator():
    '''
    Init it with a filename. File format: has N+1 lines => first line is delT, and next N lines are data point values { array[i] = value at delT*(i-1) second }
    # array[0] = time difference delT, array[1..5000] => data is presented at delT interval        
    '''
    def __init__(self, filename):
        array_in = np.genfromtxt(filename, delimiter='')
        self.array = array_in[1:]
        self.delT = array_in[0]  #time interval data were taken
        self.numData = self.array.shape[0] # 5000
        self.maxTime = self.delT * self.numData - 1
        self.minTime = 0
 
    def get_data_at_t (self, t):
        '''
        t: time point
        '''
        delT = self.delT
        array = self.array
        maxTime = self.maxTime
        minTime = self.minTime
         
        if t>maxTime or t<minTime:
            return 0
        if t%delT == 0:
            return array[int(t/delT)]
        else:
            x = t/delT
            x1 = int(np.floor(x))
            x2 = int(x1 + 1)
	    if x2*delT > maxTime:
		return array[x1]
            y = array[x1] + (array[x2] - array[x1])* (x-x1)/(x2-x1)
            return y
 

#  will imitate y = mean + amplitute * sin (omega * t + delta)
class FunctionValueGenerator(object):
	def __init__(self, filename):
		self.tracegenerator = TraceDataGenerator(filename)

	def get_value(self, t):
		return self.tracegenerator.get_data_at_t(t)

	def integrate(self, lower_limit, upper_limit):
		lowIndex = int(math.ceil(1.0*lower_limit/self.tracegenerator.delT))
		highIndex = int(math.ceil(1.0*upper_limit/self.tracegenerator.delT))
		accum = 0.0
		for i in range(lowIndex, highIndex):
			accum = accum + 0.5 * (self.tracegenerator.array[i]+self.tracegenerator.array[i+1])*self.tracegenerator.delT
		accum = accum + 0.5 * (self.tracegenerator.array[highIndex]+self.get_value(upper_limit)) * (upper_limit - highIndex*self.tracegenerator.delT)
		accum = accum + 0.5 * (self.tracegenerator.array[lowIndex]+self.get_value(lower_limit)) * (lowIndex*self.tracegenerator.delT - lower_limit)
		return accum
		

if __name__ == '__main__':
	f = FunctionValueGenerator("data_altitude.csv")
	print f.integrate(11, 13.3)
