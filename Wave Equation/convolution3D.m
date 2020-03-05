clear;

n = 128; x = linspace(-20, 20, n);
[x, y, z] = meshgrid(x, x, x);
r = sqrt(x.^2 + y.^2 + z.^2);
c = 10; dr = y(2) - y(1);

psi0 = exp(-r.^2 / 2) / sqrt((2 * pi)^3);
psi = zeros(n, n, n);
g = zeros(n, n, n);

m = n / 2;
x1D = x(m, :, m);
index = 1;
psi1D = zeros(n, 5);

tic
for t = [5, 10, 50, 100, 150]/100
    g = heaviside(t - r/c)/4/pi./r;
    psi = convn(psi0, g, 'same')*dr^3;
    psi1D(:, index) = psi(m, :, m);
    
    disp(['t = ', num2str(t), ' @ ', num2str(toc), ' seconds.']);
    index = index + 1;
end

figure; hold on;
plot(x1D, psi1D, 'LineWidth', 2);
plot(x1D, 1/4/pi./r(m,:,m), 'k--');
axis([-20 20 0 0.07]);
xlabel('x (au)');
ylabel('Scalar Potential, \phi \rightarrow 1 / 4\pir');
legend('t = 0.05','t = 0.1','t = 0.5','t = 1.0','t = 1.5','1/(4\pir)');