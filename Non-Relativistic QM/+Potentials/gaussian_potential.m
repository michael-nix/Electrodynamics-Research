function V = gaussian_potential(x,~)
%   GAUSSIAN_POTENTIAL  A gaussian potential centered at x = 0.
%
%   V = GAUSSIAN_POTENTIAL(x) returns the value of the potential for any or 
%       all values of x.
%
%   V = GAUSSIAN_POTENTIAL(x, t) returns the value of the potential at time t
%       for any or all values of x.  Only used to help with standardized
%       notation when called by other functions.
%
%   See also HARMONIC_POTENTIAL, ELECTRIC_POTENTIAL, ZERO_POTENTIAL,
%   MORSE_POTENTIAL, HARMONIC_GAUSSIAN_POTENTIAL.

V = -6 * exp(-x.^2/8);