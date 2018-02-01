import numpy as np
import math
import scipy.io as sio
from gm_Simulator import gm_Simulator
import random


def initParams(params):
    # Params in testingg.m
    params.slope = 0
    params.duty_col = np.linspace(0.5,1,11)
    params.duty_col = params.duty_col[9::-1]
    params.outer_col = - np.linspace(math.pi,0,11)
    params.outer_col = params.outer_col[1:]
    params.meter_col = math.pi + params.outer_col

    params.numParams = 9 # if != 9, then fct return 0
    params.Aleg = math.pi / 3
    params.a1 = 0.
    params.a2 = math.pi/3.
    params.b1 = 0.
    params.b2 = math.pi
    params.b3 = 0.

def fct_get_a3(a1, a2, b2, b3):
    a3_all = np.linspace(-math.pi/3., math.pi/3., 81).tolist()
    best_a3 = None
    best_dist = None
    for a3 in a3_all:
        cur_dist = fct_findMax_dist2Pi(a1, a2, a3, b2, b3)
        if (best_dist == None or cur_dist < best_dist):
            best_a3 = a3
            best_dist = cur_dist
    return best_a3


def fct_findMax_dist2Pi(a1, a2, a3, b2, b3):
    time = np.linspace(0, 2 * math.pi, 200)
    fct_func = a1 * np.sin(time / 2.) + a2 * np.sin(time + b2) + a3 * np.sin(time * 2. + b3)
    return np.max(np.abs(fct_func)) - math.pi/3.

def runSimulatorOri():
    # Init all the parameters
    class Struct(object): pass
    params = Struct()
    initParams(params)
    # Start to run the simulator with different duty_ind
    dutyf_save = np.zeros(10)
    for duty_ind in xrange(1, 10):
        #dutyf_save[duty_ind - 1] = params.duty_col[duty_ind]
        params.dutyf = dutyf_save[duty_ind - 1]
        params.interPhi = math.pi
        params.meterPhi = params.meter_col[params.phi_ind]
        params.outerPhi = params.outer_col[params.phi_ind]
        (angleChange, displacement) = gm_Simulator(params)
        print "Angle Change: " + str(angleChange) + " pi"
        print "Displacement: " + str(displacement) + "\n"

    #sio.savemat('duty_factor.mat', {'dutyf': dutyf_save})

def runSimulator(dutyFactor, interPhi, meterPhi, outerPhi, c1 = 1., a2 = 0., c2 = 1., b2=0., c3 = 1., c4 = 1., printResult = False):
    # Init all the parameters
    class Struct(object): pass
    params = Struct()
    params.numParams = 9 # if != 9, then fct return 0
    params.Aleg = [math.pi / 3, math.pi / 3]

    params.slope = 0

    params.dutyf = dutyFactor
    params.interPhi = interPhi
    params.meterPhi = meterPhi
    params.outerPhi = outerPhi

    params.a2 = [a2_1, a2_2]
    params.b2 = [b2_1, b2_2]

    params.b1 = [0., 0.]
    params.a1 = [0., 0.]
    params.a3 = [0., 0.]

    # turnning scaler for beta (leg angle)
    params.cFL = c1
    params.cHL = c2
    params.cHR = c3
    params.cFR = c4

    ######## CHANGE IT TO FALSE TO TURN OFF PLOTTING ########
    params.ifPlot = False


    # Start to run the simulator with different duty_ind
    (angleChange, displacement, unstablePenalty) = gm_Simulator(params)
    if (printResult):
        print "Angle Change: " + str(angleChange) + " pi"
        print "Displacement: " + str(displacement) + "cm"
        print "Unstable Penalty: " + str(unstablePenalty) + "\n"

    return (angleChange, displacement, unstablePenalty)

runSimulator(dutyFactor= 0.825596618927,
             interPhi=2.30325491293, meterPhi=3.02684850655, outerPhi=-1.12931199084,
             a2=0.640345443569, b2=2.87237958421,
             c1 = -0.661868934101, c2 = -0.966094411428, c3 = 0.677950475621, c4 = 0.253750011962,
             printResult=False)


