function V = harmonic_potential(x,~)

if isvector(x) && isrow(x)
    x = x.';
end

V = x.^2/2;