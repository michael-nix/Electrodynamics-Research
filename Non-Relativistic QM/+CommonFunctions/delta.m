function d = delta(x, dx)

if nargin < 2
    dx = abs(x(2)-x(1));
end

d = rect(x/dx);