function V = morse_potential(x, De)
%   MORSE_POTENTIAL  The Morse potential centered at x = 0.
%
%   V = MORSE_POTENTIAL(x) returns the value of the potential for any or 
%       all values of x.
%
%   V = MORSE_POTENTIAL(x, t) returns the value of the potential at time t
%       for any or all values of x.  Only used to help with standardized
%       notation when called by other functions.
%
%   See also HARMONIC_POTENTIAL, ELECTRIC_POTENTIAL, ZERO_POTENTIAL,
%   GAUSSIAN_POTENTIAL, HARMONIC_GAUSSIAN_POTENTIAL.

if nargin < 2
    De = 10;
end

a = sqrt(1/2/De);
V = De * (1-exp(a*x)).^2;