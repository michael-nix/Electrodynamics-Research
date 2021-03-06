clear;

% set up the workspace; adjust these variables to see how things change:
n = 1024;   xmin = -5; xmax = 5;   dt = 0.06;
x = linspace(xmin, xmax, n).';
dx = x(2) - x(1);

% V is the potential whose stationary states we want to investigate:
V = zeros(size(x));
a = -1i*dt/4/dx^2;

% set up the Hamiltonian / time-evolution operator and find eigenvectors:
H = [a*ones(n, 1), 1-2*a+1i*dt*V, a*ones(n, 1)];
H = spdiags(H, -1:1, n, n);
H = H \ conj(H);
[v, d] = eigs(H);

% figure out the energy levels, clean things up, then plot:
[ev_calc, idx] = sort(diag(abs(atan2(imag(d), real(d))/dt)));
ev_theory = (1:6).^2.' * pi^2/2 / (xmax-xmin)^2;
error = abs(ev_theory - ev_calc)./ev_theory * 100;

results = table(ev_theory, ev_calc, error);
results.Properties.VariableNames = {'Theory', 'Calculated', 'Error (%)'};

efuncs = real(v(:, idx));
for i = 1:6
    efuncs(:, i) = efuncs(:, i) / sqrt(sum(efuncs(:, i).^2*dx));
end

plot(x, efuncs, 'LineWidth', 1.5);
grid on; xlabel('x (a.u.)');
title('Calculated Eigenfunctions');
