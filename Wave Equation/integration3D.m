% setup workspace:
tmax = 150; n = 128; x = linspace(-20, 20, n);
[x, y, z] = meshgrid(x, x, x);
c = 10; dt = 0.01; dr = y(2) - y(1); s = dt^2 * c^2;

% initialize scalar potential (u) and source (f):
f = exp((-x.^2 - y.^2 - z.^2) / 2) / sqrt((2 * pi)^3);
u_now = zeros(size(x));
u_prev = zeros(size(x));
u_next = zeros(size(x));

m = n / 2;
x1D = x(m, :, m);
psi1D = zeros(n, 5);

fig1 = figure; hold on;

index = 1;
plot_times = [5, 10, 50, 100, 150];

% loop through explicit method over time of interest:
for t = 1:tmax
    u_next = 2 * u_now - u_prev + s * (6 * del2(u_now, dr, dr, dr) + f);
    u_prev = u_now;
    u_now = u_next;
    
    if any(t == plot_times)
        psi1D(:, index) = u_now(m, :, m);
        
        figure(fig1);
        plot(x1D, psi1D(:, index), 'LineWidth', 2);
        drawnow;
        
        index = index + 1;
    end
end

plot(x1D, 1/4/pi./abs(x1D), 'k--');
axis([-20 20 0 0.07]);
xlabel('x (au)');
ylabel('Scalar Potential, \phi \rightarrow 1 / 4\pir');
legend('t = 0.05','t = 0.1','t = 0.5','t = 1.0','t = 1.5','1/(4\pir)');