% setup workspace:
n = 128;   rmin = -20;   rmax = 20;

x = linspace(rmin, rmax, n);
[x, y, z] = meshgrid(x, x, x);

c = 10;   dt = 0.01;   dr = y(2) - y(1);   s = dt^2 * c^2;

% ensure stability via CFL condition:
while s / dr^2 >= 1
    dt = dt * 0.95;
    s = dt^2 * c^2;
end

% initialize scalar potential (u) and source (f):
f = exp((-x.^2 - y.^2 - z.^2) / 2) / sqrt((2 * pi)^3);
u_now = zeros(size(x));
u_prev = zeros(size(x));
u_next = zeros(size(x));

m = n / 2;
x1D = x(m, :, m);
psi1D = zeros(n, 5);

% clear up some memory and get ready to loop:
clear x y z;

plot_times = [5, 10, 50, 100, 150];
tmax = max(plot_times);
idx = 1;

fig1 = figure;   hold on;

% loop through explicit method over time of interest:
for t = 1:tmax
    u_next = 2 * u_now - u_prev + s * (6 * del2(u_now, dr, dr, dr) + f);
    u_prev = u_now;
    u_now = u_next;
    
    if any(t == plot_times)
        psi1D(:, idx) = u_now(m, :, m);
        
        figure(fig1);
        plot(x1D, psi1D(:, idx), 'LineWidth', 2);
        drawnow;
        
        idx = idx + 1;
    end
end

ncells = length(plot_times) + 1;
labels = cell(ncells, 1);
for i = 1:(ncells - 1)
    labels{i} = ['t = ', num2str(dt * plot_times(i))];
end
labels{end} = '1/(4\pir)';

plot(x1D, 1/4/pi./abs(x1D), 'k--');   axis([-20 20 0 0.07]);
xlabel('x (au)');   ylabel('Scalar Potential, \phi \rightarrow 1 / 4\pir');
legend(labels);
