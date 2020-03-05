function R = rect(x)

R = zeros(size(x));

R(abs(x) < 0.5) = 1;
R(x == 0.5) = 0.5;