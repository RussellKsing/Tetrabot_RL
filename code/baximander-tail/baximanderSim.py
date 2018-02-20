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

    params.numParams = 11 # if != 9, then fct return 0
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

def runSimulator(dutyFactor, interPhi, meterPhi, outerPhi, c1 = 1., a2_1 = 0., a2_2 = 0., a2_3 = 0., c2 = 1., b2_1 = 0., b2_2 = 0., b2_3 = 0., c3 = 1., c4 = 1., printResult = False):
    # Init all the parameters
    class Struct(object): pass
    params = Struct()
    params.numParams = 11 # if != 9, then fct return 0
    params.Aleg = math.pi / 3
    params.slope = 0

    params.dutyf = dutyFactor
    params.interPhi = interPhi
    params.meterPhi = meterPhi
    params.outerPhi = outerPhi

    params.a2 = [a2_1, a2_2, a2_3]
    params.b2 = [b2_1, b2_2, b2_3]

    params.b1 = [0., 0., 0.]
    params.a1 = [0., 0., 0.]
    params.a3 = [0., 0., 0.]

    # turnning scaler for beta (leg angle)
    params.cFL = c1
    params.cHL = c2
    params.cHR = c3
    params.cFR = c4

    ######## CHANGE IT TO FALSE TO TURN OFF PLOTTING ########
    params.ifPlot = True


    # Start to run the simulator with different duty_ind
    (angleChange, displacement, unstablePenalty) = gm_Simulator(params)
    if (printResult):
        print "Angle Change: " + str(angleChange) + " pi"
        print "Displacement: " + str(displacement) + "cm"
        print "Unstable Penalty: " + str(unstablePenalty) + "\n"

    return (angleChange, displacement, unstablePenalty)

# runSimulator(dutyFactor= 0.825596618927,
#              interPhi=2.30325491293, meterPhi=3.02684850655, outerPhi=-1.12931199084,
#              a2_1=0.640345443569, b2_1=2.87237958421, a2_2=0.78328366274, b2_2=2.81628394402,
#              c1 = -0.661868934101, c2 = -0.966094411428, c3 = 0.677950475621, c4 = 0.253750011962,
#              printResult=False)

# runSimulator(
#     dutyFactor =  0.802286750372,
#      interPhi =  2.8377008911,
#      meterPhi =  1.31511418541,
#      outerPhi =  -1.58499300215,
#      a2_1 =  0.687061588972,
#      a2_2 =  0.617446659696,
#      b2_1 =  1.54641748222,
#      b2_2 =  2.24965620214,
#      c1 =  -0.894343048431,
#      c4 =  0.0957389822278,
#      c2 =  -1.0380206765,
#      c3 =  0.341622924691,
#      printResult = False)

# runSimulator(
#     dutyFactor = 0.859385509984,
#     interPhi = 2.68174000402,
#     meterPhi = 2.15359919427,
#     outerPhi = -2.05713393347,
#     a2_1 = 0.719126912197,
#     a2_2 = 0.671344642314,
#     b2_1 = 1.79264738808,
#     b2_2 = 1.94682837911,
#     c1 = -0.877067686606,
#     c4 = -0.107887432572,
#     c2 = -1.03910473637,
#     c3 = 0.367188823342,
#     printResult = False)

# turns.. and the gait looks okay, but the robot turns back
# runSimulator(
#     dutyFactor =  0.87284995616,
#      interPhi =  2.64996121189,
#      meterPhi =  1.91151712179,
#      outerPhi =  -1.82451108116,
#      a2_1 =  0.668231307437,
#      a2_2 =  0.668075117147,
#      b2_1 =  1.68930306697,
#      b2_2 =  1.41974015163,
#      c1 =  -0.714751790867,
#      c4 =  -0.321353250046,
#      c2 =  -1.02705950861,
#      c3 =  0.684648404924,
#      printResult = False)

#manual calibration to the above gait
# does well with the sequence of actions of the legs.. 
# but need it for the body as well
# runSimulator(
#     dutyFactor =  0.87284995616,
#      interPhi =  2.64996121189,
#      meterPhi =  1.91151712179,
#      outerPhi =  -1.82451108116,
#      a2_1 =  0.668231307437,
#      a2_2 =  1.068075117147,
#      b2_1 =  1.68930306697,
#      b2_2 =  1.01974015163,
#      c1 =  -0.714751790867,
#      c4 =  0.721353250046,
#      c2 =  -1.02705950861,
#      c3 =  0.684648404924,
#      printResult = False)

# alpha2 = 0
runSimulator(
    dutyFactor = 0.9,
     interPhi = 4.01212602265,
     meterPhi = 1.18943521037,
     outerPhi = -1.41714577328,
     a2_1 = 0.701045249749,
     a2_2 = 0.905461661386,
     a2_3 = 2.5,
     b2_1 = 2.86276125191,
     b2_2 = 0.0474152993261,
     b2_3 = 0.5,
     c1 = -0.900959302491,
     c4 = 0.790483551065,
     c2 = -0.84577315574,
     c3 = 0.606085753847,
    printResult = False
    )

# mid-way simulation episode 18 - looks good
# runSimulator(
#     dutyFactor =   0.890579599126,
#      interPhi =   4.74009290177,
#      meterPhi =   0.145939537205,
#      outerPhi =   -0.0426641691804,
#      a2_1 =   0.796399083116,
#      a2_2 =   0.848591398732,
#      a2_3 =   2.82634109627,
#      b2_1 =   4.77613175603,
#      b2_2 =   1.1049991791,
#      b2_3 =   0.518174831952,
#      c1 =   -0.814722213588,
#      c4 =   0.904086508033,
#      c2 =   -0.068065407553,
#      c3 =   0.899353933873,
#     printResult = False
#     )

# runSimulator(
#     dutyFactor =  0.806335494787,
#     interPhi =  6.97765117533,
#     meterPhi =  0.514910950469,
#     outerPhi =  2.23031254112,
#     a2_1 =  0.645769527235,
#     a2_2 =  0.720532796534,
#     b2_1 =  0.544295124164,
#     b2_2 =  0.647702510908,
#     c1 =  0.437647897833,
#     c4 =  -0.757089642563,
#     c2 =  -0.813178558149,
#     c3 =  0.0354795048153,
#     printResult = False)

# runSimulator(dutyFactor= 0.864829056412,
# interPhi=2.41283940265,
# meterPhi=3.00274638295,
# outerPhi=-1.27321947551,
# a2_1=0.61837456212,
# b2_1=2.78274975643,
# a2_2=0.78328366274,
# b2_2=2.81628394402,
# c1 = -0.64826427123, 
# c2 = -0.97284655817, 
# c3 = 0.68037266484, 
# c4 = 0.23746461829,
# printResult=False)
# runSimulator(dutyFactor= 0.864829056412,
#              interPhi=2.41283940265, meterPhi=3.00274638295, outerPhi=-1.27321947551,
#              a2=0.61837456212, b2=2.88274975643,
#              c1 = -0.64826427123, c2 = -0.97284655817, c3 = 0.68037266484, c4 = 0.23746461829,
#              printResult=False)

# runSimulator(dutyFactor= 0.791440428291,
#              interPhi=3.10431216758, meterPhi=2.39017708665, outerPhi=-1.07166770839,
#              a2=0.75072157564, b2=2.36998083133,
#              c1 = -0.597754147288, c2 = -1.0, c3 = 0.226046985075, c4 = 0.472105143335,
#              printResult=False)

# * Duty Factor:  0.791440428291
#     * interPhi:  3.10431216758
#     * meterPhi:  2.39017708665
#     * outerPhi:  -1.07166770839
#     * a2:  0.75072157564
#     * b2:  2.36998083133
#     * c_FL:  -0.597754147288
#     * c_FR:  0.472105143335
#     * c_HL:  -1.0
#     * c_HR:  0.226046985075

# Duty Factor:  
#     * interPhi:  
#     * meterPhi:  
#     * outerPhi:  
#     * a2:  
#     * b2:  
#     * c_FL:  
#     * c_FR:  
#     * c_HL:  
#     * c_HR:  