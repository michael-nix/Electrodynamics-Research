function [sigmax, sigmay] = setupPML(x, dx)
    PMLwidth = 15;   m = 4;
    
    PMLwidth = PMLwidth * dx;
    sigmax_max = log(1e4)*(m+1) / PMLwidth^(m+1);
    PMLwidth = max(max(x)) - PMLwidth;
    
    sigmax = zeros(size(x));
    idx = abs(x) > PMLwidth;
    sigmax(idx) = sigmax_max * (abs(x(idx)) - PMLwidth).^m;
    sigmay = sigmax.';
    
    sigmax = sparse(sigmax);
    sigmay = sparse(sigmay);
end
