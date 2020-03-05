function R = ramp(x)

R = zeros(size(x));

R(x >= 0) = x(x >= 0);