import numpy as np
import math
import random
import scipy.io as sio
from runSimulator import runSimulator
from gm_Simulator import f_leg_act
import sys
from pprint import pprint

def initPramsVector(numParams = 9):
    '''
    Randomly init parameters
    :param numParams
    :return: 1D tuple
    '''
    #dutyf = random.uniform(0.5, 0.99)
    dutyf = random.uniform(0.6, 0.75)
    interPhi = random.uniform(0, 2 * math.pi)
    meterPhi = random.uniform(0, 2 * math.pi)
    outerPhi = random.uniform(0, 2 * math.pi)
    if (numParams == 9):
        a1 = random.uniform(0, math.pi/3.)
        a2 = random.uniform(0, math.pi / 3.)
        b1 = random.uniform(0, 2 * math.pi)
        b2 = random.uniform(0, 2 * math.pi)
        c = random.uniform(-1., 1.)
        return (dutyf, interPhi, meterPhi, outerPhi, a1, a2, b1, b2, c)
    else:
        return (dutyf, interPhi, meterPhi, outerPhi)

def policyNearPi(oriPi, epsilon):
    '''
    Generate a tuple of randomly generated policy near ori pi
    :param oriPi: the original pi
    :param epsilon: a small number
    :return: a tupleof
    '''
    R = []
    epsilonOfR = []
    for theta in oriPi:
        eps = np.random.choice([+epsilon, 0, -epsilon])
        R += [theta + eps]
        if (eps == +epsilon):
            epsilonOfR += [1]
        elif (eps == -epsilon):
            epsilonOfR += [-1]
        else:
            epsilonOfR += [0]
    return tuple(R), tuple(epsilonOfR)

# normalize the variables wrt their range
def normalizePi(oriPi, scaler):
    # Range = [dutyf, interPhi, meterPhi, outerPhi, c1, a2, c2, b2, c3, c4]
    paramsRange = [(0.5, 0.9),
                   (0, 2 * math.pi), (0, 2 * math.pi), (0, 2 * math.pi),
                   (-1., 1.), (0, math.pi/3),
                   (-1., 1.), (0, 2 * math.pi), (-1., 1.), (-1., 1.)]
    scaledPi = []
    for i in range(len(oriPi)):
        paramRangeLength = paramsRange[i][1] - paramsRange[i][0]
        scaledPi += [((oriPi[i] - paramsRange[i][0]) / paramRangeLength) * scaler]
    return scaledPi

# denormalize the variables wrt their range
def denormalizePi(scaledPi, scaler):
    # Range = [dutyf, interPhi, meterPhi, outerPhi, a1, a2, b1, b2, b3]
    paramsRange = [(0.5, 0.9),
                   (0, 2 * math.pi), (0, 2 * math.pi), (0, 2 * math.pi),
                   (-1., 1.), (0, math.pi / 3),
                   (-1., 1.), (0, 2 * math.pi), (-1., 1.), (-1., 1.)]
    oriPi = []
    for i in range(len(scaledPi)):
        paramRangeLength = paramsRange[i][1] - paramsRange[i][0]
        oriPi += [(scaledPi[i]/scaler) * paramRangeLength + paramsRange[i][0]]
    return oriPi


def evaluatePi(thePi, scaler):
    (dutyf, interPhi, meterPhi, outerPhi, c1, a2, c2, b2, c3, c4) = denormalizePi(thePi, scaler)
    # isStable = checkStable(dutyf, interPhi, meterPhi, outerPhi)
    (angleChange, displacement, unstablePercentage) = runSimulator(dutyf, interPhi, meterPhi, outerPhi,
                                                   a2=a2, b2=b2, c1=c1, c2=c2, c3=c3, c4=c4)
    w_dis = 1
    w_angle = 3
    w_stable = 7
    norm_dis = (3.5 - displacement) / 3.5
    norm_angle = angleChange
    norm_stable = 1. - unstablePercentage
    finalScore = w_dis * norm_dis + w_angle * norm_angle + w_stable * norm_stable

    return (angleChange, displacement, unstablePercentage, finalScore)



def checkStable(dutyf, interPhi, meterPhi, outerPhi):
    return True

    time = np.linspace(0, 2*math.pi, 120)
    for t in time:
        curTime = np.array([t])
        leg1 = f_leg_act(dutyf, curTime + interPhi)[0]
        leg2 = f_leg_act(dutyf, curTime + outerPhi)[0]
        leg3 = f_leg_act(dutyf, curTime)[0]
        leg4 = f_leg_act(dutyf, curTime + meterPhi)[0]
        curAct = (leg1, leg2, leg3, leg4)

        if (curAct == (1., 1., 0., 0.) or curAct == (0., 0., 1., 1.)
            or curAct == (1., 0., 1., 0.) or curAct == (0., 1., 0., 1.)
            or curAct == (1., 0., 0., 0.) or curAct == (0., 1., 0., 0.)
            or curAct == (0., 0., 1., 0.) or curAct == (0., 0., 0., 1.)
            or curAct == (0., 0., 0., 0.)):
            return False

    return True


def scaledA(A, stepSize):
    # return stepSize * (A/|A|), 1D list
    numpyA = np.array(A, dtype = float)
    sumA = sum(np.absolute(numpyA))
    if (sumA != 0):
        numpyA = (numpyA / sumA) * stepSize
    return numpyA.tolist()


def normA(A):
    # return (A/|A|), 1D list
    numpyA = np.array(A, dtype=float)
    sumA = sum(np.absolute(numpyA))
    if (sumA != 0):
        numpyA = (numpyA / sumA)
    return numpyA.tolist()

def momentumA(A, preA, stepSize, gama):
    # A, preA must be normed A
    numpyA = np.array(A, dtype=float)
    numpyPreA = np.array(preA, dtype=float)
    adjA = gama * numpyPreA + stepSize * numpyA
    return adjA.tolist()

def updatePi(thePi, adjustmentVector):
    # return pi + adjustmentVector
    return np.add(thePi, adjustmentVector).tolist()


def wrapUpParams(thePi):
    # Range = [dutyf, interPhi, meterPhi, outerPhi, a1, a2, b1, b2, b3]
    # thePhi must be Denormalized
    paramsRange = [(0.5, 0.99),
                   (0, 2 * math.pi), (0, 2 * math.pi), (0, 2 * math.pi),
                   (0, math.pi / 3), (0, math.pi / 3),
                   (0, 2 * math.pi), (0, 2 * math.pi), (0, 2 * math.pi)]

    for i in range(len(thePi)):
        # For interPhi, meterPhi, outerPhi and b1, b2, b3, if out of range
        # then wrap them up
        if (i != 0 and i != 4 and i != 5):
            upperBound = paramsRange[i][1]
            lowerBound = paramsRange[i][0]
            cycle = upperBound - lowerBound
            if (thePi[i] < lowerBound):
                thePi[i] += cycle
            elif (thePi[i] > upperBound):
                thePi[i] -= cycle

    return thePi

def setParamInBound(thePi, scaler):
    # Range = [0, scaler]
    # thePi must be normalized
    newPi = list(thePi)
    for i in range(len(newPi)):
        if (i == 0 or i == 4 or i == 5 or i == 6 or i == 8 or i == 9):
            if (newPi[i] < 0):
                newPi[i] = 0
            elif (newPi[i] > scaler):
                newPi[i] = scaler
    return tuple(newPi)



def policyGradientTrainOneEpisode(thePi, t, epsilon, object = 'Displacement', stepSize = 2., paramScaler = 7):

    print 'Step Size: ', stepSize

    # 1. Init policy
    #    initPi = [dutyf, interPhi, meterPhi, outerPhi, a1, a2, b1, b2, b3]
    thePi = normalizePi(thePi, paramScaler)
    evlIdex = 0 if (object == 'AngleChange') else 1

    # 2. Create t randomlty generated policies {R1, R2, ... Rt} near pi
    #    dicR_eval = {(Ri): (angleChange, displacement)}
    #    dicEpsilon_R = {(epsilonOfR_i): [R3, R5, ...]}
    dicR_eval = {}
    dicEpsilon_R = {}
    for i in range(t):
        R_i, epsilonOfR_i = policyNearPi(thePi, epsilon)
        curResult = evaluatePi(R_i, paramScaler)

        # Avoid unstable gait
        '''
        while (curResult[-1] < 0.0):
            R_i, epsilonOfR_i = policyNearPi(thePi, epsilon)
            curResult = evaluatePi(R_i, paramScaler)
        '''


        dicR_eval[R_i] = curResult

        dicEpsilon_R[epsilonOfR_i] = dicEpsilon_R.get(epsilonOfR_i, []) + [R_i]

    # Get the adjustment vector A for the parameters
    A = []
    for theta_ind in range(len(thePi)):
        # 3. Groping each R_i into one of three sets for each dimension n:
        S_positive = []
        S_zero = []
        S_negative = []

        for epsilonOfRs in dicEpsilon_R:
            # Ri that have a positive perturbation in dimension n
            #   * If the nth parameter of Ri is theta_n + epsilon
            if (epsilonOfRs[theta_ind] > 0):
                for R_i in dicEpsilon_R[epsilonOfRs]:
                    evaluations_Ri = dicR_eval[R_i]
                    S_positive += [evaluations_Ri[-1]]

            # Ri that have a zero perturbation in dimension n
            #   * If the nth parameter of Ri is theta_n + 0
            elif (epsilonOfRs[theta_ind] == 0):
                for R_i in dicEpsilon_R[epsilonOfRs]:
                    evaluations_Ri = dicR_eval[R_i]
                    S_zero += [evaluations_Ri[-1]]

            # Ri that have a negative perturbation in dimension n
            #   * If the nth parameter of Ri is theta_n - epsilon
            elif (epsilonOfRs[theta_ind] < 0):
                for R_i in dicEpsilon_R[epsilonOfRs]:
                    evaluations_Ri = dicR_eval[R_i]
                    S_negative += [evaluations_Ri[-1]]

        # 4. Compute the average score of the three sets
        AvgOfTheta_pos = 0
        AvgOfTheta_zero = 0
        AvgOfTheta_neg = 0
        if (len(S_positive) != 0):
            AvgOfTheta_pos =  sum(S_positive) / len(S_positive)
        if (len(S_zero) != 0):
            AvgOfTheta_zero =  sum(S_zero) / len(S_zero)
        if (len(S_negative) != 0):
            AvgOfTheta_neg =  sum(S_negative) / len(S_negative)

        # 5. Generate the adjustment for the nth parameter theta_ind
        An = AvgOfTheta_pos - AvgOfTheta_neg
        if (AvgOfTheta_zero > AvgOfTheta_pos and AvgOfTheta_zero > AvgOfTheta_neg):
            An = 0

        A += [An]

    # 6. Normalize A and  multiply it by a scalar stepSize
    #    So our adjustment will remain a fixed size each iteration
    normedA = normA(A)

    # warp up the parameters
    # denormNewPi = wrapUpParams(denormNewPi)

    return normedA

def momentumUpdatePi(A, preA, oldPi, paramScaler, stepSize, gama):
    normPi = normalizePi(oldPi, paramScaler)
    newA = momentumA(A, preA, stepSize, gama)
    newPi = updatePi(normPi, newA)
    newPi = setParamInBound(newPi, paramScaler)
    denormNewPi = denormalizePi(newPi, paramScaler)
    return denormNewPi

def decayStepsize(stepsize, iter, shrinkFactor):
    curSize =  shrinkFactor * stepsize * math.pow(1.03, -iter)
    # print 'Step Size (iter %d): %f'%(iter, curSize)
    return curSize

def getShrinkFactor(shrinkFactor, iter):
    if (iter < 3):
        return shrinkFactor
    else:
        return math.pow(0.8, iter - 2) * shrinkFactor

def getGama(initGama, iter, totalEpisodes):
    return initGama + (iter/totalEpisodes) * (0.99 - initGama)


def policyGradientTrain(t, episodes, epsilon, selfInit = False, object = 'Displacement', stepSize = 2., paramScaler = 7):
    while (True):
        thePi = initPramsVector(numParams=9)
        if (selfInit):
            thePi = selfInitParams()
        (dutyf, interPhi, meterPhi, outerPhi, c1, a2, c2, b2, c3, c4) = thePi
        if (checkStable(dutyf, interPhi, meterPhi, outerPhi)):
            break

    evlIdex = 0 if (object == 'AngleChange') else 1
    print '************* Initial ****************'
    
    # pprint (thePi)
    # pprint (normalizePi(thePi, paramScaler))
    # sys.exit()
    
    evalResults = evaluatePi(normalizePi(thePi, paramScaler), paramScaler)
    printEvaluation(thePi, evalResults, paramScaler)
    sys.exit()

    savePi = [thePi]
    saveResults = [evalResults]

    preA = np.zeros([len(thePi)], dtype=float).tolist()
    for episode in range(episodes):
        print '************* %d Episode ****************'%(episode+1)
        iter = 0
        shrinkFactor = getShrinkFactor(1.0, iter)
        curStepSize = decayStepsize(stepSize, episode, shrinkFactor)
        curBestStepSize = curStepSize

        ## Update pi
        A = policyGradientTrainOneEpisode(thePi, t, epsilon, object=object, stepSize=curStepSize,
                                         paramScaler=paramScaler)
        gama = getGama(0.5, episode, episodes)
        curBestPhi = momentumUpdatePi(A, preA, thePi, paramScaler, stepSize, gama)

        ## Evaluate
        curBestEvalResults = evaluatePi(normalizePi(curBestPhi, paramScaler), paramScaler)

        while (curBestEvalResults[-1] <= evalResults[-1] and iter < 10):
            iter += 1
            shrinkFactor = getShrinkFactor(1.0, iter)
            newStepSize = decayStepsize(stepSize, episode, shrinkFactor)

            ## Update pi
            A = policyGradientTrainOneEpisode(thePi, t, epsilon, object=object, stepSize=newStepSize,
                                                   paramScaler=paramScaler)
            gama = getGama(0.5, episode, episodes)
            newPhi = momentumUpdatePi(A, preA, thePi, paramScaler, stepSize, gama)

            newEvalResults = evaluatePi(normalizePi(newPhi, paramScaler), paramScaler)
            if (newEvalResults[-1] > curBestEvalResults[-1]):
                curBestPhi = newPhi
                curBestEvalResults = newEvalResults
                curBestStepSize = newStepSize

        # Update loop values
        thePi = curBestPhi
        evalResults = curBestEvalResults
        preA = A
        print 'this Step Size: ', curBestStepSize

        printEvaluation(thePi, evalResults, paramScaler)
        savePi += [thePi]
        saveResults += [evalResults]

        if (checkStoppingCreteria(saveResults, obj = 'Displacement')):
            break

    sio.savemat('RL_results.mat', {'parameters': savePi, 'results': saveResults})


def checkStoppingCreteria(savedResults, obj = 'Displacement'):
    if (len(savedResults) <= 5):
        return False
    if (5 < len(savedResults)):
        sumDis = 0
        for i in range(5):
            curResult = savedResults[-i-1]
            preResult = savedResults[-i-2]
            sumDis += abs(curResult[-1] - preResult[-1])
            print'sumDis: '+ str(sumDis)
        if (len(savedResults) >= 30 and sumDis/5.0 < 0.005):
            return True
        elif (len(savedResults) < 30 and sumDis == 0.):
            return True
        else:
            return False

    return False


def printEvaluation(thePi, evalResults, paramScaler):
    (angleChange, displacement, unstablePercentage, penaltiedScore) = evalResults
    (dutyf, interPhi, meterPhi, outerPhi, c1, a2, c2, b2, c3, c4) = thePi
    print '------------------------------------------------------'
    print 'Parameters: '
    print '    * Duty Factor: ', dutyf
    print '    * interPhi: ', interPhi
    print '    * meterPhi: ', meterPhi
    print '    * outerPhi: ', outerPhi
    #print '    * a1: ', a1
    print '    * a2: ', a2
    #print '    * b1: ',  b1
    print '    * b2: ', b2
    print '    * c_FL: ', c1
    print '    * c_FR: ', c4
    print '    * c_HL: ', c2
    print '    * c_HR: ', c3

    print
    print "Angle Change: " + str(angleChange) + " pi"
    print "Displacement: " + str(displacement) + " cm"
    print "Unstable Percentage: " + str(unstablePercentage * 100) + '%'
    print "Final Score: " + str(penaltiedScore) + '\n'

def selfInitParams():
    dutyf =  random.uniform(0.70, 0.95)
    legPhaseShift = 30./100

    interPhi = 3.40803280698#math.pi
    meterPhi = 1.96412224952#legPhaseShift * 2 * math.pi
    outerPhi = -1.48728484168#math.pi + legPhaseShift * 2 * math.pi

    '''
    a1 = 0.162100870114
    a2 = 0.162100870114
    b1 = 5.72129100099
    b2 = 5.69980412613
    b3 = 4.36691547831
    '''
    a2 = 0.825796303018#random.uniform(0, math.pi / 3.)
    b2 = 2.48508330853#random.uniform(0, 2 * math.pi)

    c1 = -random.uniform(-0.2, 1) # cFL
    c2 = -random.uniform(-0.2, 1) # cHL
    c3 = random.uniform(-0.2, 1) # cHR
    c4 = random.uniform(-0.2, 1) # cFR


    '''
    dutyf =  0.6
    interPhi = 3.4
    meterPhi = 3.4
    outerPhi = 0.1
    a1 = 0.803112592501
    a2 = 1.01
    b1 = 5.72129100099
    b2 = 3.2
    b3 = 0
    '''

    '''
    * Duty Factor:  0.5
    * interPhi:  pi
    * meterPhi:  pi
    * outerPhi:  0
    * a2:  pi/3
    * b2:  pi
    * b3:  0
    '''

    return (dutyf, interPhi, meterPhi, outerPhi, c1, a2, c2, b2, c3, c4)

policyGradientTrain(t=30, episodes=80, epsilon=0.001, selfInit = True, object = 'AngleChange', stepSize = 0.5, paramScaler = 7)


# thePi = [0.5, 2.4, 3.8, 1, 0.5, math.pi/3, 2.1 + 2 * math.pi, 3.15, 0.17]
# print policyGradientTrainOneEpisode(thePi, 1, 0.0001, object = 'Displacement', stepSize = 2, paramScaler = 7)





