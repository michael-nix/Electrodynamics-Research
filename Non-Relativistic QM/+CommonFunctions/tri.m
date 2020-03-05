function T = tri(x)

T = zeros(size(x));

index = abs(x) < 1;
T(index) = 1-abs(x(index));