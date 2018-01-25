## Proposal

- previously, handcrafted a sidewinding gait
- but it can't use the environment intelligently

### What was the gait?
- 2 wave system phased properly
	- several modulations of phase and amplitude contribute to "good templates"

- global behaviors are represented using low dimensional parameterized expressions

- systematic discovery of templates include the following steps from geometric mechanics:
	- generate the height function
	- draw curve over the zero set

### Questions proposed to be answered
- find templates in heterogenous environment
- find templates for legged robots with cyclic gaits
- analysis and visualization on sidewinding snakes, 4-legged salamanders and 6-legged force feedback robots

### Maths

- Original equation
$$
	\xi = A(\alpha)\dot{\alpha} + \Gamma(\alpha)\dot{p}
$$
- If inertial effects are negligeble, 
$$
	\xi = A(\alpha)\dot{\alpha}
$$

- we represent the shape of a complex locomotor as linear combination of only two shape basis functions
$$
	\alpha = w_1\beta_1 + w_2\beta_2
$$
- the kinetic reconstruction equation becomes
$$
	\xi = A(w)\frac{\partial \alpha}{\partial w}\dot{w}
$$

### Symmetric Gaits

1. lateral limbs are always out of phase to each other
2. phase difference between FL and HL = phase difference between FR and HR

- 2 parameters
	- duty factor: % of stride for each leg to be on the ground
	- leg phase shift: % of stride, fore footfall follows the hind on the same side
- Hildebrand diagram
	- diagonal sequence
	- lateral sequence

### Technical Work

- identify optimal templates for snakes, tetrapods, hexapods
- extend geometric-mechanics to apply systems with intermittent contacts
- identify low dimensional representation for motion in heterogeneous environment
- hypothesis: modification to these templates will support heterogeneous environment
- Bio experiments to determine how snakes, walking salamanders modulate their templates in heterogeneous env

#### Hybrid Contacts

- Intermittent contacts
- activation variable: $\delta \in [0,1]$ 0: no contact, 1: full contact
- A, now a function of $w$ and $\delta$
- Also, $\delta$ now a function of $\alpha$

#### Basis Functions

- reduce dimension by using shape basis functions
- modulation upon the shape basis function can improve locomotor performance in heterogeneous env
- computing shape basis functions:
	1. shape basis optimization algorithm - iteratively improves the choice of shape basis function
	2. Learn from biological systems
	3. learn empirically from robot experiments

- combine basis function design with hybrid contacts

#### Machine Learning

- explore motion beyond the templates
- templates make assumptions
	- relaxing these assumptions, gives modes of motion
	- explored modes of motion using policy gradient RL
