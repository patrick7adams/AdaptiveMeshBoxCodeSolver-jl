# Error Analysis

Testing a few knobs to turn to get a sense of where the error lies in the current code. Sections will be given by the parameter that I change.

## Default values

These default values will be used for all parameters, other than the one being changed:

| Parameter | Value |
| --------- | ----- |
| $\delta$  | $0.1$ |
| Curvature | $200$ |
| $\sigma$  | $0.05$|
| $x_0$       | $0.85$|

## Parameter descriptions

$\delta$ is the relative size of the boundary strip against the distance between boundary points.

Curvature controls the number of discretization points along the boundary.

$\sigma$ is the strength of the gaussian function that we are using.

$x_0$ is the x-position of the gaussian.

## Method description

We will be fitting the mesh to the following function:
$$\Delta u (x, y) = \left( \frac{(x - x_0)^2 + (y - y_0)^2}{\sigma^4} - \frac{2}{\sigma^2} \right) u(x, y)$$
Where $u(x, y)$ is given by:
$$ u(x, y) = e^{-\frac{(x-x_0)^2 + (y-y_0)^2}{2\sigma^2}}$$

With default $y_0 = 0$.

The mesh that we are fitting has a boundary given by the unit circle. All domain quadratures have quadrature order $4$, and all boundary quadratures have quadrature order $6$. 

Error will be given via the max error and the L2 error for each variation.

## Default values

With the default values, the errors are given below:

$$L_2 \text{ error} = 8.57*10^{-8}$$
$$\text{Max error} = 2.41*10^{-6}$$

## Mesh Separation Distance

This parameter describes the relative threshold (against quad size) of using the singular quad integration method as opposed to the nonsingular quadrature method. Default value is $1.5$. 

| Threshold | $L_2 \text{ error}$ | Max error |
| --------- | ------------------ | --------- |
| $2.0$     | $7.37*10^{-8}$     | $2.41*10^{-6}$ |
| $1.9$     | $7.33*10^{-8}$     | $2.39*10^{-6}$ |
| $1.8$     | $7.33*10^{-8}$     | $2.39*10^{-6}$ |
| $1.7$     | $7.34*10^{-8}$     | $2.39*10^{-6}$ |
| $1.6$     | $7.50*10^{-8}$     | $2.39*10^{-6}$ |
| $1.5$     | $8.57*10^{-8}$     | $2.41*10^{-6}$ |
| $1.4$     | $1.11*10^{-7}$     | $2.43*10^{-6}$ |
| $1.3$     | $1.19*10^{-7}$     | $2.43*10^{-6}$ |
| $1.2$     | $1.20*10^{-7}$     | $2.43*10^{-6}$ |
| $1.1$     | $1.29*10^{-7}$     | $2.43*10^{-6}$ |
| $1.0$     | $8.01*10^{-7}$     | $2.03*10^{-6}$ |


## Delta

Variance in this type of error will show the accuracy of the triangle potential calculations as opposed to the quad potential calculations.

| $\delta$ | $L_2 \text{ error}$ | Max error |
| --------- | ------------------ | --------- |
| $0.0$    | $7.48*10^{-8}$     | $1.21*10^{-6}$ |
| $0.05$    | $7.46*10^{-8}$     | $2.01*10^{-6}$ |
| $0.1$     | $8.57*10^{-8}$     | $2.41*10^{-6}$ |
| $0.2$     | $8.76*10^{-8}$     | $2.51*10^{-6}$ |
| $0.4$     | $7.73*10^{-8}$     | $7.92*10^{-7}$ |
| $0.8$     | $1.06*10^{-7}$     | $2.81*10^{-6}$ |
| $1.6$     | $1.14*10^{-7}$     | $3.12*10^{-6}$ |
| $3.2$     | $2.52*10^{-7}$     | $2.26*10^{-5}$ |
| $6.4$     | $2.45*10^{-7}$     | $4.24*10^{-6}$ |
| $12.8$     | $9.74*10^{-7}$     | $2.80*10^{-5}$ |

## Curvature

| Curvature | $L_2 \text{ error}$ | Max error |
| --------- | ------------------ | --------- |
| $50$    | $5.53*10^{-7}$     | $9.72*10^{-6}$ |
| $100$    | $1.15*10^{-7}$     | $2.36*10^{-6}$ |
| $200$     | $8.57*10^{-8}$     | $2.41*10^{-6}$ |
| $400$     | $6.99*10^{-8}$     | $1.65*10^{-6}$ |

## $\sigma$

| $\sigma$ | $L_2 \text{ error}$ | Max error |
| --------- | ------------------ | --------- |
| $0.05$    | $8.57*10^{-8}$     | $2.41*10^{-6}$ |
| $0.1$    | $1.86*10^{-7}$     | $3.88*10^{-6}$ |
| $0.2$     | $1.31*10^{-7}$     | $5.46*10^{-6}$ |
| $0.4$     | $3.21*10^{-7}$     | $2.17*10^{-6}$ |
| $0.8$     | $3.86*10^{-8}$     | $4.11*10^{-7}$ |

## $x_0$

| $x_0$ | $L_2 \text{ error}$ | Max error |
| --------- | ------------------ | --------- |
| $0.95$    | $6.83*10^{-7}$     | $3.13*10^{-5}$ |
| $0.85$    | $8.57*10^{-8}$     | $2.41*10^{-6}$ |
| $0.7$    | $7.70*10^{-8}$     | $7.57*10^{-7}$ |
| $0.5$     | $7.53*10^{-8}$     | $7.45*10^{-7}$ |
| $0.0$     | $1.14*10^{-7}$     | $8.19*10^{-7}$ |

## Note on the Rc vs. c coeffs

After the fix was applied (multiplying by R to get coeffs), accuracies improved very slightly. This is expected, as the coeffs in the other basis were still able to indicate where refinement was necessary; they just didn't indicate it properly. Generally, all errors here can be considered slightly reduced, but the same patterns have seemingly arisen.

I'd like to run another test of this, but I'd first like to work on the necessary speedups; namely, fixing separateMesh is my top priority