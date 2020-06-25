clear;

% setup workspace:
n = 128;   rmin = -20;   rmax = 20;

x = linspace(rmin, rmax, n);
[x, y, z] = meshgrid(x, x, x);
r = sqrt(x.^2 + y.^2 + z.^2);

c = 10;   dr = y(2) - y(1);

% initialize scalar potential (psi), green's function (g), and source (f):
f = exp(-r.^2 / 2) / sqrt((2 * pi)^3);
psi = zeros(n, n, n);
g = zeros(n, n, n);

m = n / 2;
x1D = x(m, :, m);
psi1D = zeros(n, 5);

% clear up some memory and get ready to loop:
clear x y z;

plot_times = [5, 10, 50, 100, 150] / 100;
idx = 1;

tic;
% loop through convolutions over time of interest:
for t = plot_times
    g = heaviside(t - r/c) / 4/pi./r;
    psi = convn(f, g, 'same') * dr^3;
    psi1D(:, idx) = psi(m, :, m);
    
    disp(['t = ', num2str(t), ' @ ', num2str(toc), ' seconds.']);
    idx = idx + 1;
end

% plot all of the results:
ncells = length(plot_times) + 1;
labels = cell(ncells, 1);
for i = 1:(ncells - 1)
    labels{i} = ['t = ', num2str(plot_times(i))];
end
labels{end} = '1/(4\pir)';

figure;    hold on;

plot(x1D, psi1D, 'LineWidth', 2);   plot(x1D, 1/4/pi./r(m,:,m), 'k--');
axis([rmin rmax 0 0.07]);

xlabel('x (au)');   ylabel('Scalar Potential, \phi \rightarrow 1 / 4\pir');
legend(labels);
