## Kinematic Gait synthesis for snake robots

### Aim of the paper

- to analyze snake kinematics
- shape trajectory data is high dimensional
	- hence apply modal decomposition to reduce dimension
- PCA, etc simplifies instantaneous shapes
	- over time analysis is a challenge
- not many have studied 3D simulation
	- divide into horizontal/vertical - 2 sets of shape trajectory
	- conditional basis array factorization
		- high dim array to low dim parametrized representation
	- CDAF: good for over time analysis

### Biological Snakes
- turning behavior
	1. Differential turning: gradual turn + forward motion
	2. Reversal turning: change in orientation + small turn radius
- recording snake movement details

### Undulations in orthogonal planes
- recent: vertical undulations are independent of horizontal undulations
- N $=$ 8 (joint angles), T $=$ 128 (frames)
- N x T tables for horizontal and vertical planes

Identify low dim parametrized representation of high order array
- parametrize shape trajectory of individual behavior + parametrize change in body undulations across different behavior

- N: spatial direction (angles)
- T: temporal direction (# frames)
- G: behavioral direction (# sequences)

### CBAF
- builds on array factorization
- modeling snake locomotion as the linear combination of a few spatiotemporal modes

- model decomposition technique: PCA - requires first order array
- multidim structure of high order data array
	1. HOSVD
	2. ALS
	Both computationally expensive

### Factorization with conditional bases

- Condition bases: chosen based on domain knowledge or candidate bases
- Using condition bases -
	1. Basis fitting problem to basis selection problem
		- choosing a subset of columns