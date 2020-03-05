function u = initeven(x, x0, p, sig)
%   INITEVEN  Initializes a gaussian wavefunction symmetric about x0.
%
%   u = INITEVEN(x, x0, p) creates a gaussian wavefunction centered at x0
%       with initial momentum given by p.
%
%   See also INITODD.

if nargin == 1
    x0 = 0;
end

if nargin < 4
    sig = 1;
end

if isvector(x)
    n = 1;
else
    n = ndims(x);
end

u = exp(-(x - x0).^2/2/sig^2) / (sqrt(2*pi)*sig)^n;

if nargin >= 3
    u = u .* exp(1i*p*x);
end