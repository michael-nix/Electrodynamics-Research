function H = heaviside(x)

H = zeros(size(x));

H(x > 0) = 1;
H(x == 0) = 0.5;