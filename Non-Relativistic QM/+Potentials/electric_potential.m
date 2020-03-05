function V = electric_potential(x,~)
%   ELECTRIC_POTENTIAL  The electric potential centered at x = 0.
%
%   V = ELECTRIC_POTENTIAL(x) returns the value of the potential for any or
%       all values of x.
%
%   V = ELECTRIC_POTENTIAL(x, t) returns the value of the potential at time
%       t for any or all values of x.  Only used to help with standardized
%       notation when called by other functions.
%
%   Absorbing potentials are added at either end of the domain of the
%   ProblemSpace to assist with simulation.
%
%   See also HARMONIC_POTENTIAL, GAUSSIAN_POTENTIAL, ZERO_POTENTIAL,
%   MORSE_POTENTIAL, HARMONIC_GAUSSIAN_POTENTIAL.

V = -1./abs(x);

if min(x) <= -50 && max(x) >= 50
    PMLwidth = round(length(x) * 0.05);
    dx = x(2) - x(1);
    PMLwidth = ceil(PMLwidth * dx);
    sigma_max = log(1e8)*5 / PMLwidth^5;
    PMLwidth = max(x) - PMLwidth;
    
    idx = abs(x) > PMLwidth;
    V(idx) = V(idx) - 1i*sigma_max*(abs(x(idx)) - PMLwidth).^4;
end