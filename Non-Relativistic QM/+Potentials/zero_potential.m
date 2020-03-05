function V = zero_potential(x,~)
%   ZERO_POTENTIAL  A potential that is zero within the given ProblemSpace.
%   If used in conjunction with the ThomasAlgorithm solver it can behave as
%   an infinite square potential well.  This is true if it's also used with
%   any solver with Dirichlet boundary conditions such that the function is
%   zero at the boundary.
%
%   V = ZERO_POTENTIAL(x) returns the value of the potential for any or all
%       values of x.
%
%   V = ZERO_POTENTIAL(x, t) returns the value of the potential at time t
%       for any or all values of x.  Only used to help with standardized
%       notation when called by other functions.
%
%   See also HARMONIC_POTENTIAL, ELECTRIC_POTENTIAL, MORSE_POTENTIAL,
%   GAUSSIAN_POTENTIAL, HARMONIC_GAUSSIAN_POTENTIAL.

V = zeros(size(x));