import numpy as np
from math import sin, cos, pi
from pprint import pprint
import sys

def computeBodyVelocity(alpha, d_alpha, beta, d_beta, hind, L, Lleg, activation, K, F):
    '''
    * parameters are all numpy arrays
    :param alpha:   1x3 1D numpy vector, body angle
    :param d_alpha: 1x3 1D numpy vector, derivative of alpha
    :param beta:    1x4 1D numpy vector, leg angle
    :param d_beta:  1x4 1D numpy vector, derivative of beta
    :param hind:    a constant (hind = 3.0)
    :param L:       a constant (L = 5.5)
    :param Lleg:    a constant (L = 7.0)
    :param activation:  1x8 1D numpy vector, activation coefficients, 1.2 or 0
                        [bodyx4, front left, front right, back left, back right]
    :param K:       constant numpy array 3x3 2D list
    :param F:       constant 1x3 1D numpy vector, F = [0,0,0]

    :return xi:
    '''

    # compute the body velocity of the head frame
    # compute all the frames relative to the head frame
    n = alpha.shape[0] # n = 3

    # Get all the modules frame relative to head
    (g, gFL, gFR, gHL, gHR) = framesInHead(alpha, beta, hind, L, Lleg)

    ########################
    # Compute omega1
    ########################

    # activation: [body, front left, front right, back left, back right]
    # g is a [4 x 3 x 3] 3D numpy array
    # omega1 = activation(1) * (myAdjoint(inv(g{1}))).'*K*activation(1);
    omega1 = activation[0] * np.dot(np.transpose((myAdjoint(np.linalg.inv(g[0])))), activation[0] * K)

    for i in xrange(1, n): # i = 1,2
        # omega1 = omega1+activation(i)*(myAdjoint(inv(g{i}))).'*K*myAdjoint(inv(g{i}));
        omega1 += omega1Helper(activation[i], g[i], K)


    # Front left
    '''
    -- Matlab code --
    omega1 = omega1 + activation(n+1) * (myAdjoint(inv(gFL))).'*K*myAdjoint(inv(gFL));
    omega1 = omega1 + activation(n+2) * (myAdjoint(inv(gFR))).'*K*myAdjoint(inv(gFR));
    omega1 = omega1 + activation(n+3) * (myAdjoint(inv(gHL))).'*K*myAdjoint(inv(gHL));
    omega1 = omega1 + activation(n+4) * (myAdjoint(inv(gHR))).'*K*myAdjoint(inv(gHR));'''

    omega1 += omega1Helper(activation[n+1], gFL, K)
    omega1 += omega1Helper(activation[n+2], gFR, K)
    omega1 += omega1Helper(activation[n+3], gHL, K)
    omega1 += omega1Helper(activation[n+4], gHR, K)

    ########################
    # Compute omega2
    ########################

    # Compute spatial Jacobian, front legs, hind legs
    (J1, Jf, Jh) = spatialJacobian(alpha, hind, L)
    omega2 = np.zeros([3, 1])

    for i in xrange(1, n):
        # for i = 2:n
        # omega2 += activation(i)*(myAdjoint(inv(g{i}))).'*K*myAdjoint(inv(g{i}))*[J1(:,1:i-1),zeros(3,n-i)]*d_alpha;
        padJ = np.concatenate((J1[:, 0:i], np.zeros([3,n-i])), axis=1)
        omega2 += omega2Helper(activation[i], g[i], K, padJ, d_alpha)

    # Legs
    omega2 += omega2Helper(activation[n+1], gFL, K, Jf, d_beta[0])

    omega2 += omega2Helper(activation[n+2], gFR, K, Jf, d_beta[1])

    padJ = np.concatenate((J1[:,0:hind - 1], np.zeros([3, n - hind + 1]), Jh), axis = 1)
    d_betaA = np.concatenate((d_alpha,[d_beta[2]]))
    omega2 += omega2Helper(activation[n+3], gHL, K, padJ, d_betaA)

    d_betaA = np.concatenate((d_alpha,[d_beta[3]]))
    omega2 += omega2Helper(activation[n+4], gHR, K, padJ, d_betaA)
    omega2 += np.reshape(F, [3,1])

    ########################
    # xi
    #
    xi = np.dot(-np.linalg.inv(omega1), omega2)
    return xi





def omega1Helper(activation, g, K):
    # return: activation * (myAdjoint(inv(g1))).'*K*myAdjoint(inv(g2));
    return np.dot(activation * np.dot(np.transpose((myAdjoint(np.linalg.inv(g)))), K),
                         myAdjoint(np.linalg.inv(g)))

def omega2Helper(activation, g, K, padJ, d_alpha):
    # return: activation * (myAdjoint(inv(g1))).'*K*myAdjoint(inv(g2)) * [J,zeros(3,n-i)]*d_alpha;
    firstPart = omega1Helper(activation, g, K)
    return np.reshape(np.dot(np.dot(firstPart, padJ), d_alpha),[3,1])


def myAdjoint(g):
    # compute adjoint opteration from g
    row1 = np.concatenate((g[0, 0:2], [g[1, 2]]), axis=0)
    row2 = np.concatenate((g[1, 0:2], [-g[0, 2]]), axis=0)
    row3 = np.array([0, 0, 1])
    ad = np.array((row1, row2, row3))
    return ad

# input parameters:
# alpha: amplitude
def framesInHead(alpha, beta, hind, L, Lleg):
    n = alpha.shape[0]
    g = np.zeros((4,3,3)) # 4x3x3 np array
    g[0] = np.eye(3,3) # module 0 configuration: identity in itself
    F = np.array([[1,0,L],[0,1,0],[0,0,1]]) # move forward
    # pprint (g)
    for i in xrange(1, n-1):
        # new module frame
        # g{i}=g{i-1}*F*R(alpha(i-1))*F
        # pprint (np.dot(g[i-1], F))
        # pprint (rotationFrame(alpha[i-1]))
        # pprint (alpha[i-1])
        g[i] = np.dot(np.dot(np.dot(g[i-1], F), rotationFrame(alpha[i-1])), F)
    # pprint (g)
    # sys.exit()
    F = np.array([[1,0,L/2.],[0,1,0],[0,0,1]])
    # g{n-1}=g{n-2}*F*R(alpha(n-2))*F*F;
    g[n-1] = np.dot(np.dot(np.dot(np.dot(g[n-2], F), rotationFrame(alpha[n - 1])), F), F)

    F = np.array([[1, 0, L / 3.], [0, 1, 0], [0, 0, 1]])
    # g{n}=g{n-1}*F*F;
    g[n] = np.dot(np.dot(g[n-1], F),F)

    F = np.array([[1, 0, L], [0, 1, 0], [0, 0, 1]], dtype=np.float32)
    # gFL = g{1}/F*R(beta(1))*[eye(2),[Lleg;0];0 0 1];
    gFL = np.dot(np.dot(np.dot(g[0], np.linalg.inv(F)), rotationFrame(beta[0])), np.array([[1,0,Lleg],[0,1,0],[0,0,1]]))
    # pprint (gFL)
    # sys.exit()
    # gFR = g{1}/F*R(beta(2))*[eye(2),[Lleg;0];0 0 1];
    gFR = np.dot(np.dot(np.dot(g[0], np.linalg.inv(F)), rotationFrame(beta[1])), np.array([[1,0,Lleg],[0,1,0],[0,0,1]]))
    # gHL = g{hind-1}*F*R(beta(3))*[eye(2),[Lleg;0];0 0 1];
    gHL = np.dot(np.dot(np.dot(g[hind - 2], F), rotationFrame(beta[2])), np.array([[1,0,Lleg],[0,1,0],[0,0,1]]))
    # gHR = g{hind-1}*F*R(beta(4))*[eye(2),[Lleg;0];0 0 1];
    gHR = np.dot(np.dot(np.dot(g[hind - 2], F), rotationFrame(beta[3])), np.array([[1, 0, Lleg], [0, 1, 0], [0, 0, 1]]))

    return (g, gFL, gFR, gHL, gHR)


def spatialJacobian(alpha, hind, L):
    '''
    Spatial Jacobian computes the spatial manipulator Jacobian with the head
    Module defined as the spatial reference frame

    :param alpha: 1 x n, joint angles
    :param hind: constant
    :param L: constant
    :return: J1, Jf, Jh: 3 x n
    '''

    # compute joint frame position:
    # note the first frame is the head not the joint
    n = alpha.shape[0] + 1 # 3
    q = np.zeros([2, n - 1]) # joint position
    g_joint = jointsInHead(alpha, L)

    for i in xrange(0, n-1):
        q[:, i] = g_joint[i+1, 0:2, 2]

    # Get all the joint positions, and we now the rotational axis is [0 0 1]
    # Now, we construct the spatial jacobian

    q = np.concatenate(([q[1,:]], [-q[0,:]])) # Position part in the spatial Jacobian
    J1 = np.concatenate((q, np.ones([1, n-1])))

    # Front legs
    q0 = g_joint[0, 0:2, 2]
    q0 = np.array([[q0[1]],[-q0[0]]])

    Jf = np.concatenate((q0, [[1]]))

    qHind = g_joint[hind - 1, 0:2, 2]
    qHind = np.array([[qHind[1]],[-qHind[0]]])
    Jh = np.concatenate((qHind, [[1]]))

    return (J1, Jf, Jh)


def jointsInHead(alpha, L):
    # compute all the joint frame, defined as the proximal end of the next
    # Module at the joint
    n = alpha.shape[0] + 1 # number of modules, 4
    g = np.zeros((n + 1, 3, 3)) # g is [5,3,3]
    g[0] = np.eye(3,3) # Put the head frame at the identity
    F = np.array([[1,0,L],[0,1,0],[0,0,1]]) # Move forward half link length
    g[0] = np.dot(g[0], np.linalg.inv(F)) # The virtual 0th joint
    # Put 0 joint back

    for i in xrange(1, n-1):
        g[i] = np.dot(np.dot(np.dot(g[i-1], F), F), rotationFrame(alpha[i-1]))

    g[n-1] = np.dot(np.dot(g[n-2], F), rotationFrame(alpha[n-2]))
    F = np.array([[1,0,L/3.],[0,1,0],[0,0,1]])
    g[n] = np.dot(np.dot(g[n-1], F), F)

    return g


def rotationFrame(a):
    return np.array([[cos(a), -sin(a), 0],
                     [sin(a),  cos(a), 0],
                     [     0,       0, 1]])


def test():
    alpha = np.array([-0.2052,0,0])
    d_alpha = np.array([0.8959,0,0])
    beta = np.array([1.5708, -1.2489, 2.5506, -2.2512])
    d_beta = np.array([3.8859, 0.5753, 1.3045, 0.4611])
    hind = 3
    L = 5.5
    Lleg = 7
    K = np.array([[1.,0.,0.],
                  [0.,2.,0.],
                  [0.,0.,1.]])
    activation = [1.2,1.2,0,0,0,1,1,1]
    F = np.array([0.,0.,0.])
    #xi = computeBodyVelocity(alpha, d_alpha, beta, d_beta, hind, L, Lleg, activation, K, F)
    (g, gFL, gFR, gHL, gHR) = framesInHead(alpha, beta, hind, L, Lleg)
    print g
    print
    print gFL
    print
    print gFR
    print
    print gHL
    print
    print gHR



#test()

