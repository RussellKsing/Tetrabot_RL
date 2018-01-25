## Policy Gradient Reinforcement Learning for Fast Quadrapedal Locomotion

### What are they doing?

- Optimizing quadrapedal trot gait
- Policy gradient RL to automatically search the set of possible parameter -> to find the dfastest walk.

>what is gait?

### Robot description

- Sony Aibo Robot
- Robocup
- Hand tuning parameterized gait
- Need ML to automate the search for good parameters
- Need ML to dfinding optimal control policies

### What is the norm?

- Gait is determined by a series of joint positions for the three joints in each of its legs
- These days, we deal with trajectory of Aibo's four feet in 3D space and then perform inverse kinematics that converts the trajectories to joint angles
- We choose shape of the loci of the feet and the parameters of these loci
	- trapezoid loci
	- eleptical loci

### What is a Gait Trot

- Each pair of diagonally opposite legs in phase with each other and perfectly out of phase with the other two

### 12 Parameters

1. front locus (3 param: x, y, $\theta$)
2. rear locus (3 param)
3. locus length
4. locus skew multiplier
5. height of the front of the body
6. height of the rear of the body
7. time each foot takes to move through the locus
8. fraction of time each foot spends on the ground

- speed is the sole objective function
- considering each possible set of parameter assignment as defining open-loop policy that can be executed by the robot
> what is an open-loop policy?
- we estimate policy's gradient and then follow it towards a local minima. Difficulty:
	- no functional form
	- sampling is expensive

### Efficient method to estimate the policy gradient

- Q-learning not possible because:
	- no notion of state
	- open-loop control
	- no MDPs

The paper's approach:
1. Parameter vector: $\pi = ${$\theta_1, \ldots, \theta_N$} (N=12, $\theta_i$ are the parameters)
2. We need to estimate $\frac{\partial Obj}{\partial\theta_i}$. How?
	- generating random policies {$ R_1, \ldots, R_t $}
	- each $R_i = ${$ \theta_1+\Delta_1,\ldots,\theta_N+\Delta_N $ }
	  where $\Delta_i = +\epsilon, 0, -\epsilon$ relative to $\theta_i$
3. Evaluate speed at each $R_i$ -> which they call the score
4. Estimate $\frac{\partial Obj}{\partial\theta_i}$ by dividing into three groups
	- g1: all $R_i$'s with parameter 1 to be $\theta_i - \epsilon_i$
	- g2: all $R_i$'s with parameter 1 to be $\theta_i$
	- g3: all $R_i$'s with parameter 1 to be $\theta_i + \epsilon_i$
5. Compute average _score_ of these groups and call them $Avg_{-\epsilon,n}, Avg_{0,n}, Avg_{+\epsilon,n},$
6. These averages gives us the benifit of altering the $n^{th}$ parameter by $-\epsilon, 0, +\epsilon$. See figure 3 in paper for better understanding.
7. We use these scores to construct an adjustment vector A of size N. where
$$
	A_n = 0, \text{ if } Avg_{0, n} > Avg_{+\epsilon, n} \text{ and } Avg_{0, n} > Avg_{-\epsilon, n}
$$
$$	
	A_n = Avg_{+\epsilon, n} - Avg_{-\epsilon, n}, \text{ otherwise}
$$
8. Normalize A, multiply by step size $\eta$
9. Add A to $\pi$

### Results

- Velocity sometimes decreases with increasing interations but overall they achieve the highest velocity as compared to previous gaits. 
- These decrease in velocity can be attributed to
	- noise
	- large step size
- Reason for success are
	- their policy gradient algo is particularly effective in learning this task
	- the gait implementation itself is superior to the previous approaches