function u = initodd(x, x0, p, sig)
%   INITODD  Initializes a gaussian wavefunction anti-symmetric about x0.
%
%   u = INITODD(x, x0, p) creates a gaussian wavefunction centered at x0
%       with initial momentum given by p.
%
%   See also INITEVEN.

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

u = (x - x0).*exp(-(x - x0).^2/2/sig^2) / (sqrt(2*pi)*sig)^n;

if nargin >= 3
    u = u .* exp(1i*p*x);
end