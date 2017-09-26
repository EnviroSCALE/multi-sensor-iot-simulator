import heapq;
from pdb import set_trace as BP
import random;
import Queue;
import json;
from Function import FunctionValueGenerator;
import math;
import numpy as np;
from scipy.optimize import minimize

from time import sleep
#BP = lambda: None
############################################################################
# Read from config #
############################################################################
try:
    with open("sensor_config.json", 'r') as f:
        c = json.load(f)
except IOError:
    print IOError
    print("Error reading from config file: using default configuration")


class LogPrint:
	def __init__(self, obj, time):
		self.obj = obj
		self.time = time

	def __str__(self):
		return "%f:: %s" % (self.time, self.obj)



class Reading:
    def __init__(self, sensor_name, sensing_time, size, reading_value):
        self.sensor_name = sensor_name
        self.sensing_time = sensing_time
        self.size = size
        self.reading_value = reading_value

    def __repr__(self):
        return 'Sensor::%s, Time:%f, Value:%f' % (self.sensor_name, self.sensing_time, self.reading_value)



# MRH : sensor ja ja data collect korse, r trace e ja ja data ase, duita analyze kore error mapa lagbe. sheta ekhane hoy.
class DataAnalyzer(object):
    def __init__(self, f, q, sensor):
        self.f = f
        self.q = q
        self.sensor = sensor
        self.t = []
        self.v = []

    def sort_q_into_array(self):
        # MRH : kon time e kon reading, sheta niye dui array te rakhlam
        while not self.q.empty():
            reading = self.q.get()
            self.t.append(reading.sensing_time)
            self.v.append(reading.reading_value)

    # MRH : ekhane besh joghonno code krsi, sorry. -_-
    def get_error(self):
        self.sort_q_into_array()
        min_of_max_t = min(self.t[len(self.t)-1], self.sensor.f.tracegenerator.maxTime)
        time = self.t[0]
        deltaTime = 0.1
        accum = 0.0
        atIndex = 0
        while time < min_of_max_t:
            if self.t[atIndex] == self.t[atIndex+1] or not (time >= self.t[atIndex] and time <= self.t[atIndex+1]):
                atIndex = atIndex + 1
                if atIndex == len(self.t):
                    break
                continue
            val1 = 1.0 * (self.v[atIndex+1]-self.v[atIndex]) * (time - self.t[atIndex]) / (self.t[atIndex+1] - self.t[atIndex]) + self.v[atIndex]
            val2 = self.sensor.f.get_value(time)
            accum = accum + math.fabs(val1 - val2)
            time = time + deltaTime
            if time > self.t[atIndex + 1]:
                atIndex = atIndex + 1
                if atIndex == len(self.t):
                    break
        return accum

class Callable:
    def call(sim):
        raise Exception('Call of %s is not implemented' % self)


class Sensor(Callable):

    def __init__(self, name, readlatency, period, size, weight, filename, sf):
        self.name = name
        self.readlatency = readlatency
        self.period = period
        self.size = size
        self.weight = weight
        self.flag = False
        self.f = FunctionValueGenerator(filename)
        q = Queue.Queue()
        self.analyzer = DataAnalyzer(self.f, q, self)
        self.accumulated_reading = 0.0
        self.accumulated_squared_reading = 0.0
        self.total_reading_count = 0
        self.accumulated = 0.0
        self.accumulated_squared = 0.0
        self.accumulated_diff = 0.0
        self.max_in_interval = -99999.0
        self.min_in_interval = +99999.0
        self.prev_value = -1.0
        self.sf = sf

    def __repr__(self):
        return '%f \t Sensor::%s' % (sim.simclock, self.name)

    def clear(self):
        self.accumulated = 0.0
        self.accumulated_squared = 0.0
        self.total_reading_count = 0
        self.max_in_interval = -99999.0
        self.min_in_interval = +99999.0

    def set_period(self, period):
        if period == 0:
            return
        self.period = period

    def call(self, sim):
        # MRH : sensor er read latency ase. tai ekbar event e khali kisu kore na, read start hoy.
        # MRH : arek bar reading neya sesh hoy. note the add_event times. ekhane shomossa hote pare.
        if self.flag:
            self.flag = False
            #print 'Time %f sensor %s reading completed' % (sim.simclock, self)
            sim.add_event(sim.simclock + self.period - self.readlatency, self)
        else:
            self.flag = True
            sim.read_queue = sim.read_queue + self.size

            # MRH : ei value ta pabe koi theke? trace theke... sheta FunctionValueGenerator (f) handle kore.

            value = self.f.get_value(sim.simclock)
            reading = Reading(self.name, sim.simclock, self.size, value)
            sim.total_generated_data += self.size
            sim.readings_queue.put(reading)
            #print 'Time %f reading sensor %s current queue %d %f' % (sim.simclock, self, sim.read_queue, value)
            sim.add_event(sim.simclock + self.readlatency, self)

            # MRH : er porer gula holo so far seen data theke linearly variance calculate korar jonno stats gather
            alpha = self.sf
            self.accumulated_reading = self.accumulated_reading * (1-alpha) + value * (alpha)
            self.accumulated_squared_reading = self.accumulated_squared_reading * (1-alpha) + value * value * alpha
            self.accumulated = self.accumulated + value
            self.accumulated_squared = self.accumulated_squared + value * value
            self.total_reading_count = self.total_reading_count + 1
            if value > self.max_in_interval:
                self.max_in_interval = value
            if value < self.min_in_interval:
                self.min_in_interval = value
            if self.prev_value > 0:
                self.accumulated_diff = self.accumulated_diff * (1-alpha) + math.fabs(value - self.prev_value) * alpha
            self.prev_value = value

class Uploader(Callable):

    # here, bandwidth, up_time and down_time ... these are parameters that will vary from
    # network to network, place to place

    def __init__(self, period, bandwidth, upload_rate, up_time, down_time):
        self.period = period
        self.bandwidth = bandwidth
        self.upload_rate = upload_rate
        self.last_uploadtime = 0
        self.last_uploaded = 0

        self.up_time = up_time
        self.down_time = down_time

        self.failed = False
        self.flag = False
        self.currently_uploading = Queue.Queue()

    def __repr__(self):
        return 'Uploader'

    def call(self, sim):
        if self.failed:
                # MRH : jodi network failure hoy, taileo to queue ta exhaust korte hobe. ideally ei popped data niye kisu ekta korar kotha... sheta korbe FailureHandler. ei call() fn theke data gula nibe, niye kisu ekta korbe. amra apatoto FailureHandler e eshob kisu rakhi nai. pathaite na parle baad. locally o store kortesi na.
            while not self.currently_uploading.empty():
                reading = self.currently_uploading.get()

            return

        if self.flag:
            # MRH : ekhane if-else er duita block one at a time kaaj kore. like sensor, upload er khetreo... upload shuru korlam, sathe sathe to pouche jabe na. majhe kisu time lagbe. kottuk lagbe sheta network dependent. bandwidth er upor depend kore. ekhane ekta block upload initiate kora simulate kore, arekta upload sesh howa.
            self.flag = False
            # MRH : upload sesh. tai currently_uploading theke ber kore feltesi.
            while not self.currently_uploading.empty():
                # MRH : ber kore, ei j reading from a certain sensor ashlo, sheta original trace k koto valo approximate korlo shetar jonno analysis lagbe, tai ei reading take j sensor theke ashche tar analyzer e rekhe dilam.
                reading = self.currently_uploading.get()
                
                #print reading
                for s in sim.sensors:
                    if s.name == reading.sensor_name:
                        s.analyzer.q.put(reading)
                        break
                #print '----- %f delay encountered' % (sim.simclock - reading.sensing_time)
                sim.total_delay = sim.total_delay + (sim.simclock - reading.sensing_time)
                sim.total_sent = sim.total_sent + 1

            sim.add_event(sim.simclock + self.period, self)
        else:
            bytes_to_upload = max(self.upload_rate * sim.simclock - self.last_uploaded, 0)
            bytes_to_upload = math.floor(1.0*(bytes_to_upload - c["params"]["alpha"]) / c["params"]["beta"])
            bytes_to_be_uploaded = 0

            while (bytes_to_be_uploaded < bytes_to_upload) and not sim.readings_queue.empty():
                reading = sim.readings_queue.get()
                self.currently_uploading.put(reading)
                bytes_to_be_uploaded = bytes_to_be_uploaded + reading.size + c["overhead"]
            bytes_to_upload = bytes_to_be_uploaded
            # MRH : koto byte upload kora jabe sheta determine kore oi onujayi internal queue (readings_queue) theke niye uploading queue te niye rekhe dilam. hishab ta ektu dekhe bujhte hbe.

            # MRH : kotokhon lagbe upload korte, bw dekhe ber korlam
            upload_duration = 1.0 * bytes_to_upload / self.bandwidth
            sim.read_queue = max(sim.read_queue - bytes_to_upload, 0)
            # MRH : total koto bytes ami pathaisi sheta update korlam
            self.last_uploaded = self.last_uploaded + bytes_to_upload

            #print 'Time %f UPLOADING %d bytes, remaining %d in queue' % (sim.simclock, bytes_to_upload, sim.read_queue)
            #print 'So far uploaded %d bytes' % self.last_uploaded
            self.flag = True
            sim.add_event(sim.simclock + upload_duration, self)



class PeriodUpdater(Callable):
    def __init__(self, sensors, time_gap, uploader, choice):
        self.sensors = sensors
        self.interval = time_gap
        self.uploader = uploader
        self.choice = choice

    def __repr__(self):
        return "Period updater: %s" % (self.choice)

    def update(self, sim, choice):
        # MRH : nanan choice er jonno nanan code.
        # MRH : konta kon scheme sheta eqn dekhe bujha lagbe. but maxm of these amra baad diya disi, khali adaptv tai valo lagse amader kase r ki. still, inspect korte chaile kora jete pare.
        
        
        
        sim.add_event(sim.simclock + self.interval, self)
        
        # Equal weights approach
        if choice == 1:
            k = len(sim.sensors)            
            alpha = c["params"]["alpha"]
            beta = c["params"]["beta"]
            rate = self.uploader.upload_rate
            T = self.uploader.period
            if rate*T == alpha:
                return
            for i in range(0, k):
				#ART
                pi = max(0, 1.0 * k * (sim.sensors[i].size + c["overhead"]) * beta * T / (rate * T - alpha))
                sim.sensors[i].set_period(pi)
                
        # Accumulated CV approach        
        if choice == 2:
            k = len(sim.sensors)
            total = 0.0
            w = []
            for i in range(0, k):
                e_x2 = sim.sensors[i].accumulated_squared_reading
                e_x = sim.sensors[i].accumulated_reading
                var = math.sqrt(math.fabs(e_x2 - e_x*e_x))
                cv = var / e_x
                w.append(cv)
                total = total + cv
            alpha = c["params"]["alpha"]
            beta = c["params"]["beta"]
            rate = self.uploader.upload_rate
            T = self.uploader.period
            if rate*T == alpha:
                return
            if total == 0:
                return
            for i in range(0, k):
                wt = w[i] / total
                if wt == 0:
                    return
                pi = max(0, 1.0 * (sim.sensors[i].size + c["overhead"]) * beta * T / ( (rate * T - alpha) * wt ))
                sim.sensors[i].set_period(pi)
        
        # Separated CV approach      
        if choice == 3:
            k = len(sim.sensors)
            total = 0.0
            w = []
            for i in range(0, k):
                if sim.sensors[i].total_reading_count == 0:
                    return
                e_x2 = sim.sensors[i].accumulated_squared / sim.sensors[i].total_reading_count
                e_x = sim.sensors[i].accumulated / sim.sensors[i].total_reading_count
                #ART
                if e_x == 0:
                    return
                var = math.sqrt(math.fabs(e_x2 - e_x*e_x))
                cv = var / e_x
                w.append(cv)
                total = total + cv
            alpha = c["params"]["alpha"]
            beta = c["params"]["beta"]
            rate = self.uploader.upload_rate
            T = self.uploader.period
            if rate*T == alpha:
                return
            if total == 0:
                return
            for i in range(0, k):
                wt = w[i] / total
                if wt == 0:
                    return
                pi = 1.0 * (sim.sensors[i].size + c["overhead"]) * beta * T / ( (rate * T - alpha) * wt )
                #ART
                pi = max(0, pi)

                sim.sensors[i].set_period(pi)
            for s in sim.sensors:
                s.clear()
                
         # Total Error approach              
        if choice == 4:
            k = len(sim.sensors)
            total = 0.0
            w = []
            for i in range(0, k):
                frac = sim.sensors[i].accumulated_diff/(sim.sensors[i].f.integrate(sim.simclock-self.interval, sim.simclock))
                w.append(frac)
                total = total + frac
            alpha = c["params"]["alpha"]
            beta = c["params"]["beta"]
            rate = self.uploader.upload_rate
            T = self.uploader.period
            if rate*T == alpha:
                return
            if total == 0:
                return
            for i in range(0, k):
                wt = w[i] / total
                if wt == 0:
                    return
                pi = 1.0 * (sim.sensors[i].size +  c["overhead"]) * beta * T / ( (rate * T - alpha) * wt )
                sim.sensors[i].set_period(pi)
            for s in sim.sensors:
                s.clear()     
                   
        # loss function = ln (f)        
        if choice == 5:
			k = len(sim.sensors)            
			alpha = c["params"]["alpha"]
			beta = c["params"]["beta"]
			rate = self.uploader.upload_rate
			T = self.uploader.period
			if rate*T == alpha:
				return
			for i in range(0, k):
				#ART
				pi = max(0, 1.0 * (1/sim.sensors[i].weight) * (sim.sensors[i].size + c["overhead"]) * beta * T / (rate * T - alpha))
				sim.sensors[i].set_period(pi)        
		
		# loss function = square  		
        if choice == 6:
			k = len(sim.sensors)
			
			denom = 0.0     
			for i in range(0, k):
				denom += math.sqrt(sim.sensors[i].weight*sim.sensors[i].size)
				  
			alpha = c["params"]["alpha"]
			beta = c["params"]["beta"]
			rate = self.uploader.upload_rate
			T = self.uploader.period
			if rate*T == alpha:
				return
			for i in range(0, k):
				#ART
				pi = max(0, 1.0 * (denom/math.sqrt(sim.sensors[i].weight)) * (math.sqrt(sim.sensors[i].size) + c["overhead"]) * beta * T / (rate * T - alpha))
				sim.sensors[i].set_period(pi)        
                
        # No update			
        if choice == 7:
			k = len(sim.sensors)
			for i in range(0, k):
				#ART
				pi = 0
				sim.sensors[i].set_period(pi)           

    def fn(self, x):
        alpha = c["params"]["alpha"]
        beta = c["params"]["beta"]
        rate = self.uploader.upload_rate
        T = self.uploader.period
        A = (rate*T - alpha)/(beta*T)
        tot = 0.0
        for s in sim.sensors:
            tot = tot + (s.size +  c["overhead"]) / s.period
        d = tot - A
        return 0.5 * d * d

    def call(self, sim):
        self.update(sim, self.choice)

class RateUpdater(Callable):
    def __init__(self, time_gap, uploader):
        self.interval = time_gap
        self.uploader = uploader

    def __repr__(self):
        return "Rate updater"

    def call(self, sim):
        M = c["params"]["M"]
        t = sim.simclock
        if (M == t):
            return
        rate = self.uploader.upload_rate
        u = self.uploader.last_uploaded
        #assert  (c["params"]["D"] - u) * 1.0 / (M - t)) > 0
        if ((c["params"]["D"] - u) * 1.0 / (M - t) > 0):
			self.uploader.upload_rate = min((c["params"]["D"] - u) * 1.0 / (M - t), c["params"]["max_rate"])
        #self.uploader.upload_rate()
        #insufficient budget
        #assert  self.uploader.upload_rate < 0,  "Negative upload rate"
        # MRH : max rate k exceed kora jabe na. naile dekha jabe second e 50000 byte pathano jay, sheta korar jonno sampling freq almost infinity nite hobe... sheta korle to baash
        sim.add_event(sim.simclock + self.interval, self)

class FailureHandler(Callable):
    def __init__(self, uploader):
        self.uploader = uploader
        self.flag = False

    # MRH : fail korle ki korbo sheta ekhane lekha thakbe. ekhane code e only event add kora hoise.
    # MRH : ashole kemne failed data handle korbo (like storing locally or something...) sheta r kisu lekha nai
    # MRH : network ki failed naki na, sheta simulate kora hoitese 'flag' variable diye
    # MRH : eta uploader.failed = True / False diye uploader k janano hoitese. shei onujayi uploader either pathabe or pathabe na.
    def call(self, sim):
        if self.flag:
            duration = random.expovariate(1.0 / self.uploader.up_time)
            self.uploader.failed = False
            self.flag = False
            sim.add_event(sim.simclock, self.uploader)
            sim.add_event(sim.simclock + duration, self)
        else:
            duration = random.expovariate(1.0 / self.uploader.down_time)
            duration = 0
            # MRH : ekhane duration = 0 na dile crash kore. keno kore Allah janen.
            self.uploader.failed = True
            self.flag = True
            sim.add_event(sim.simclock + duration, self)

    def __str__(self):
		return "FailureHandler Class: Time is %f" % (sim.simclock)
        
        
        
class Simulator:
    def __init__(self, seed, choice, sf, update_rate):
        self.simclock = 0.0
        self.total_generated_data = 0.0
        self.event_queue = []
        self.readings_queue = Queue.Queue()
        self.read_queue = 0
        self.total_delay = 0.0
        self.total_sent = 0
        self.sf = sf
        # MRH : choice holo kon scheme use korbe sheta
        self.choice = choice
        self.update_rate = update_rate
        random.seed(seed)

    def set_endtime(self, time):
        self.endtime = time

    def init_scene(self):
        # MRH : sensor gulo nilam
        self.sensors = []
        num_sensors = len(c["sensors"])
        for i in range(0, num_sensors):
            s1 = Sensor(c["sensors"][i]["name"], c["sensors"][i]["readlatency"], c["sensors"][i]["period"],
                    c["sensors"][i]["size"], c["sensors"][i]["gamma"], c["sensors"][i]["datafile"], self.sf)
            self.sensors.append(s1)

        self.bought_data = c["params"]["D"]
        rate = 1.0 * self.bought_data / end_time
        
        
        #assert rate > 1
        #        
        
        # MRH : joto data kinsi vaag total time = rate. ei rate bojay rakhar try korbo shob shomoy
		
        upload_interval = c["upload"]["period"]
        u = Uploader(upload_interval, c["upload"]["bw"], rate, c["upload"]["up_time"], c["upload"]["down_time"])
        self.u = u
        f = FailureHandler(u)
        # MRH : data fail hoile ki korbe

        period_update_interval = c["update_interval"]["period"]
        p = PeriodUpdater(self.sensors, period_update_interval, u, self.choice)
        # MRH : eto second por por periodically sampling frequency update hobe. etar moddhe amader scheme gula code kora ase. based on choice, kono ekta scheme follow korbe

        rate_update_interval = c["update_interval"]["rate"]
        r = RateUpdater(rate_update_interval, u)
        # MRH : rateupdater holo rate update kore, mane target rate. mane ami 20 bytes/sec pathanor try kortesi always, majhe jdi kisukkhon off thake amar system, network down thakar jonno, taile to ekhon ami 22 bytes/sec pathate parbo (say). eita kore rate updater

        # shob event add kortesi
        for s in self.sensors:
            self.add_event(0, s)
        self.add_event(0, u)
        self.add_event(0, f)
        self.add_event(c["update_interval"]["period"], p) #100
        # MRH : rate update = True hoilei khali rate updater k add korbo
        if self.update_rate:
            self.add_event(c["update_interval"]["rate"], r)

    def add_event(self, time, event):
        ##print "Event generated: " + str(event) + " "  + str(time) 
        assert time >= 0
        heapq.heappush(self.event_queue, (time, event))

    def run(self):
        ###
        while len(self.event_queue) > 0:
            time, event = heapq.heappop(self.event_queue)
            
            #print time
            #print event
            if time > self.endtime:
            # MRH : sesh, shob stats hishab kortese
                self.tot = 0.0
                k = len(self.sensors)
                for i in range (0, k):
                    if i == k - 1:
                        self.prominent_error = self.sensors[i].analyzer.get_error()
                        self.tot = self.tot + self.prominent_error
                    else:
                        self.tot = self.tot + self.sensors[i].analyzer.get_error()
                self.avg_delay = 1.0 * self.total_delay / self.total_sent
                #print self.u.last_uploaded
                self.utilization = 100.0 * self.u.last_uploaded / c["params"]["D"]
                break

            self.simclock = time
            event.call(self)

    def get_data(self):
        # MRH : shobi hishab kora ase, khali ret kortese
        return self.tot, self.prominent_error, self.avg_delay, self.utilization
import csv   

if __name__ == '__main__':

    # MRH : bought data
    #D = 2000
    #D = 250000
    
    #c["params"]["D"] = D

# MRH 123 - seed
# 1 - choice, 1 means adaptive sampling, 2,3,4 means other schemes. shutro dekhe bujhte hobe
# 0.8 is the smoothing factor.
# True - whether or not I will use the 'update rate' equation
    for D in range(100000, 2000000, 100000): 
        c["params"]["D"] = D
        results = []
        sim = Simulator(c["seed"], c["choice"], c["sf"], c["is_update"])
        end_time = c["params"]["M"]
        sim.set_endtime(end_time)
        sim.init_scene()
        sim.run()
#self.tot, self.prominent_error, self.avg_delay, self.utilization, choice, time (s), data (bytes)
#0          1                       2               3               4           5       6
        result = sim.get_data() 
        result = list(result)
        result.append(c["choice"])
        result.append(c["params"]["M"])
        result.append(c["params"]["D"])
        #self.tot, self.prominent_error, self.avg_delay, self.utilization
        print 'Tot error %f' % result[0]
        print 'Error %f' % result[1]
        print 'Average delay %f' % result[2]
        print 'Utilization %f' % result[3]
        print sim.total_generated_data
        print c["params"]["D"]
        with open(c["filename"], 'a') as f:
            writer = csv.writer(f)
            writer.writerow(result)
