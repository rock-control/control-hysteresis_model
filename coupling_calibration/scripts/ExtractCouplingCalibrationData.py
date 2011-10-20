import os
import sys
import numpy as np
import pylab as pl
import math

sys.path.append(os.path.dirname(os.path.expanduser('../')))

import signalsmooth
import string

MOTOR = 0
g_angle = 0 

SMOOTH = 300
TIME = 0
wPOS = 1
wVEL = 2
wACC = 3
mPOS = 4
mVEL = 5
dPOS = 6
dVEL = 7
mCUR = 8

def differentiate(x, t, smooth_factor):
    n = np.shape(x)[0]
    xdot = np.zeros(n)
    for i in xrange(1,n):
	xdot[i] = (x[i]-x[i-1])/(t[i]-t[i-1])
    xdot = signalsmooth.smooth(np.array(xdot),smooth_factor,'bartlett')
    return xdot

def DEG(x):
    return x * 180.0 / math.pi

def RAD(x):
    return x * math.pi / 180.0

def getData(FILE):
    data = np.loadtxt(FILE, skiprows=1)
    time      = (data[:,1] - data[0,1])*0.001;
    g_angle = data[0,MOTOR*4 + 4]

    delta_pos = data[:,MOTOR*4 + 4] - data[:,MOTOR*4 + 3]
    wheel_pos = data[:,MOTOR*4 + 4] - g_angle
    motor_pos = data[:,MOTOR*4 + 3] - g_angle
    motor_cur = data[:,MOTOR*4 + 2]

    wheel_vel = differentiate(wheel_pos, time, SMOOTH)
    motor_vel = differentiate(motor_pos, time, SMOOTH)
    wheel_acc = differentiate(wheel_vel, time, SMOOTH)

    delta_vel = differentiate(delta_pos, time, SMOOTH)

    return np.transpose(np.vstack((time, wheel_pos, wheel_vel, wheel_acc, motor_pos, motor_vel, delta_pos, delta_vel, motor_cur)))

def plot(data):
    pl.figure()
    pl.plot(data[:,TIME], data[:,wPOS ], ',', label='Wheel Position (rad)')
#    pl.plot(data[:,TIME], data[:,mPOS ], ',', label='Motor Position (rad)')
#    pl.plot(data[:,TIME], data[:,mCUR ] / 1000.0, ',', label='Motor Current (A)')
    pl.plot(data[:,TIME], DEG(data[:,dPOS]), ',', label='Deflection (deg)')
    pl.plot(data[:,TIME], DEG(data[:,dVEL]), ',', label='Def Vel (deg)')
    pl.legend()
    pl.grid()

def computeTorque(data, s):
    M_tip = 0.981 + 1.099 + 0.148
    M_rod = 0.115
    D_tip = 0.2
    D_rod = 0.1

    mg_Sin_theta = -(D_tip * (M_tip + M_rod) * 9.81) * np.sin(data[:,wPOS] ) 
    ma = - (M_tip * D_tip**2 +  M_rod *  D_rod**2) * data[:,wACC] 
    torque = ma + mg_Sin_theta

    f_torque = (-torque[s[2]] + torque[s[4]])/2.0

    torque[    :s[1]]   -= f_torque
    torque[s[1]:s[3]]   += f_torque
    torque[s[3]:    ]   -= f_torque

    return torque 

def writeToFile(FILE, data, torque):
    fout=open(FILE, 'w')
    for i in xrange(0, np.shape(data)[0]):
	    fout.write(str(data[i,TIME]));
	    fout.write(' ');
	    fout.write(str(torque[i]));
	    fout.write(' ');
	    fout.write(str(data[i,dPOS]));
	    fout.write(' ');
	    fout.write(str(data[i,dVEL]));
	    fout.write('\n');
    fout.close()

def findStationaryParts(data, min_length = 1000, max_length = 5000):
    indexes = np.zeros((0,2)) 
    i = 0
    while i < np.shape(data)[0] - min_length:
    	for j in xrange(i + 1, i + max_length):
	    if data[i] != data[j] or j == np.shape(data)[0]-1:
		break
	if j - i >= min_length and j - i <= max_length:
	    indexes = np.vstack((indexes,[i,j]))
	    i = j
	i=i+1    

#    indexes[1,:] = indexes[2,:]
#    indexes[2,:] = indexes[3,:]
#    indexes[3,:] = indexes[4,:]
#    indexes[4,:] = indexes[5,:]
    return indexes

def main():
    if (len(sys.argv) != 3):
	print "Wrong number of arguments"
	print "Correct usage : python ExtractTorqueCalibrationData.py [PATH] [WHEEL NR.]"
	sys.exit()

    if not os.path.isfile(sys.argv[1]+'hbridge.txt'):
        print "hbridge.txt :File not found"
        print "Extracting log... "
    	os.system("pocolog " + sys.argv[1] + "lowlevel.0.log -s hbridge.status_motors > " + sys.argv[1] + "hbridge.txt")
    	print "Finished"
    	print "Loading data"

    global MOTOR
    MOTOR = int(sys.argv[2])
    data = getData(sys.argv[1] + 'hbridge.txt')
    plot(data)

    range_stay = findStationaryParts(data[:,wPOS], 1000, 4000)

    mid_stay = (range_stay[:,0]+range_stay[:,1]) / 2
    for i in xrange(0, np.shape(mid_stay)[0]):
	print data[mid_stay[i], TIME]

    torque = computeTorque(data, mid_stay)
#    pl.plot(data[:,TIME], torque, ',', label='Torque(Nm)')
#    pl.legend()


    data   = data  [range_stay[0,1]:,...]
    torque = torque[range_stay[0,1]:]

#    print max(data[:,dPOS])
#    print min(data[:,dPOS])
#    print ((max(data[:,dPOS]) + min(data[:,dPOS])) / 2.0)
    data[:,dPOS] -= ((max(data[:,dPOS]) + min(data[:,dPOS])) / 2.0)

    pl.figure()
    pl.plot( DEG(data[:,dPOS]), torque, ',')
    pl.grid()

    
#    pl.figure()
#    pl.plot( torque, data[:,mCUR]/1000.0, ',')
#    pl.grid()
    
    pl.show()

    writeToFile(sys.argv[1] + 'hbridge_out.txt', data, torque)

if __name__ == "__main__":
    main()
