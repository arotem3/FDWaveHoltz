## README
This repo implements the WaveHoltz algorithm using the finite difference method.
The code included can be used to solve the wave equation in the time domain:

$$u_{tt} = \Delta u + f(x, t).$$

Or solve the wave equation in the frequency domain:

$$\Delta \hat{u} + \omega^2 \hat{u} = f(x).$$

In the latter case, the WaveHoltz iteration is used:
$$(\hat{u}, \hat{v})^{(n+1)} = \Pi (\hat{u}, \hat{v})^{(n)},$$
where
$$\Pi(\hat{u}, \hat{v}) = \frac{2}{T} \int_0^T \left(\cos(\omega t) - \frac{1}{4}\right) (u, v) dt, \quad T = \frac{2\pi}{\omega},$$
$$u_t = v,$$
$$v_t = \Delta u - f(x) \cos(\omega t),$$
$$u(0) = \hat{u}, v(0) = \hat{v}.$$
