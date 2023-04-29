## Wave equation 

### Problem setup
$$\partial_{tt} u(x,t) = c^2 \partial_{xx} u(x,t)$$

$$u(x,0)= \exp{\left(-\left(\frac{x-μ_1}{σ}\right)^2\right)} + \exp{\left(-\left(\frac{x-μ_2}{σ}\right)^2\right)}$$

$$u'(x,0)=\frac{d}{dx}u(x,0)$$

$$x\in[0,\pi] \quad t\in[0,10]$$

$$\mu_1=\frac{\pi}{3} \quad \mu_1=\frac{2\pi}{3}$$

$$\sigma=0.5 \quad c=1.5$$

### Solution
 Explicit method          |  Static heatmap  
:-------------------------:|:-------------------------:
  [![name](https://github.com/dynamic-queries/LinearPDEs/blob/main/figures/wave.gif)](https://github.com/dynamic-queries/LinearPDEs/blob/main/figures/wave.gif) | [![name](https://github.com/dynamic-queries/LinearPDEs/blob/main/figures/contour_wave.png)](https://github.com/dynamic-queries/LinearPDEs/blob/main/figures/contour_wave.png)

---

## Schrodinger's equation 
With double well potential

### Problem setup
$$i\bar{h} \partial_{t} \psi(x,t) = H\psi(x,t)=\left(-\frac{\bar{h}^2}{2m}\Delta_x +V(x,t)\right)\psi(x,t)$$

$$V(x,t)= V_B \left[-\frac{1}{4}\left(\frac{z}{z_0}\right)^2+\frac{1}{64}\left(\frac{z}{z_0}\right)^4 \right] + eE_zsin(\omega t)$$

$$\psi(x,0)=\frac{\phi_1(x)+\phi_2(x)}{\sqrt{2}}$$

$$H\phi_1(x)=E_1\phi_1(x) \textrm{ and } H\phi_2(x)=E_2\phi_2(x)$$ 

$$E_1,E_2 \in \mathbb{R}$$

$$z_0=\frac{a}{4\sqrt{2}} \quad a=1$$

$$E_z=0 \quad V_B = 1 \quad \omega=200\pi\times10^{12}$$

### Solutions

This system is known to be stiff. To demonstrate this we juxtapose solutions from a 5th order adaptive-explicit vs a 2nd order implicit integrator.

Explicit method           |  Implicit method | Static  heatmap |
:-------------------------:|:-------------------------:|:-------------------------|
[![name](https://github.com/dynamic-queries/LinearPDEs/blob/main/figures/SE_explicit.gif)](https://github.com/dynamic-queries/LinearPDEs/blob/main/figures/SE_explicit.gif)  |  [![name](https://github.com/dynamic-queries/LinearPDEs/blob/main/figures/SE_implicit.gif)](https://github.com/dynamic-queries/LinearPDEs/blob/main/figures/SE_implicit.gif) | [![name](https://github.com/dynamic-queries/LinearPDEs/blob/main/figures/contour_SE.png)](https://github.com/dynamic-queries/LinearPDEs/blob/main/figures/contour_SE.png)
