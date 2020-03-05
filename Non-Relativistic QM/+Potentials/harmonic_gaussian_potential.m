function V = harmonic_gaussian_potential(x,~)

V = x.^2/2+5*exp(-5*x.^2);