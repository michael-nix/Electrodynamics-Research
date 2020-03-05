function H = heaviside(x)
if nargin < 1
    x = linspace(-pi, pi, 128);
end
H = zeros(size(x));
H(x > 0) = 1;
H(x == 0) = 0.5;