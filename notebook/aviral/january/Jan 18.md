## Simplifying Gait Design via Shape Basis optimization

### Premise

- Locomotion is complex
- Geometric mechanics (limited to 2 joints)
> what is geometric mechanics?
- use shape basis functions - intuitive
> what are shape basis functions?
- _Aim_: To find appropriate shape basis functions
- _How?_: shape basis aproximation algorithm

- They show that shapes of systems that thave many DOFs can be approzimated as a linear combination of two shape basis functions

### Geometric Mechanic

- Input: shape basis function
- Output: Design Gait for complex systems
- The next set of equations describe the mechanism of GM:
	- eqution of motion:
$$
	\xi = A(\alpha)\dot{\alpha}
$$
	- $\xi$ - body velocity
	- $\alpha$ - joint angles
	- $A(\alpha)$ - local connection: relates shape velocity ($\dot{\alpha}$) to body velocity ($\xi$)
	- shape variables: joint angles
	- group variables: body velocity

### Connection Vector Fields and Height Fucntion

- 2 joints: $\alpha \in R^2$
- $A(\alpha)$ - 3 x 2
- $\xi \in R^3$ (x, y, $\theta$)
- connection vector field is each row of the local connection matrix, a vector of dimension 2. eg: $\vec{x} = x_1 \hat{\alpha_1} + x_2 \hat{\alpha_2}$
- if $\vec{x}$ and $\vec{\alpha}$ are aligned, we get a large increase in $x$

### Gait

- it is the cyclic shape change
	- represented as closed curve in the shape space: $\partial\phi$
- displacement resulting from gait can be represented as:
$$
	\int_{\partial\phi} A(\alpha)d\alpha = \iint_{\phi}\triangledown \times A(\alpha)d\alpha_1d\alpha_2 \text{ (stokes theorem)}
$$
- $\triangledown \times A(\alpha)$ is called the height function
- details and calculation of stokes law - in notes

- Aim is to determine a gait that produces largest displacement per cycle
	- It is done by finding the zero set of the height function, which has the largest area integral 
> why does the zero set with largest area integral produce largest displacement?

### Average Body Frame

- with properly chosen body frame, called _minimum perturbation coordinates_,
$$
	\int \text{body velocity } = \text { actual displacement}
$$

- we use average body frame - (x, y, $\theta$)
	- as effective as min. pertub. coord.
	- easier to compute
- basically take average of all the link's (x, y, $\theta$)  wrt. $0^{th}$ link.

### Shape Basis Function