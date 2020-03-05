clear all; import CommonFunctions.* Fields.* Potentials.*; %#ok<CLALL>
% TESTQMSOLVER2D  Script to test QuantumSolver2D and run some experiments.
%   
%   Modify the electric potential, pot, and magnetic field, mag, to run
%   various pre-defined test simulations.
%   
%   Some pre-defined potentials are @zero_potential, @harmonic_potential,
%   @gaussian_potential, @electric_potential, @morse_potential.
%   
%   For complete list of available potentials try:
%       what Potentials
%   
%   See also QSPARAMETERS, QUANTUMSOLVER1D, QUANTUMSOLVERS.

pot = @zero_potential;

want2plot = true;

% 2^14 @ 128 = < 40 s.
% 2^14 @ 256 = < 270 s.
% 2^14 @ 512 = < 1100 s.
% 2^14 @ 1024 = < 4450 s.
tmax = 2^14; n = 256;

rmax = 50;
x = linspace(-rmax, rmax, n);
y = linspace(-rmax, rmax, n);
dx = abs(x(2) - x(1));
dy = abs(y(2) - y(1));
[x, y] = meshgrid(x, y);

% r = sqrt(x.^2+y.^2);
% theta = atan2(y, x);
sig = 1;

u0 = initeven(x, 0, 0, 2).*initeven(y, 0, 0, 2);
u0 = u0 / sqrt(sum(sum(abs(u0).^2 * dx * dy)));

% double slit experiment setup:
% u0 = initeven(x, -7).*initeven(y, 0, 0, 3).*exp(1i*x*5);
% xrange = (x > -5 & x < -3);
% pot = 1000*(1 - ((y > 0.25 & y < 1) + (y > -1 & y < -0.25))).*xrange;
% pml = zeros(size(x));   idx = abs(x) > (10 - 15*dx);
% pml(idx) = log(1e4)*5/(15*dx)^5 * (abs(x(idx)) - (10 - 15*dx)).^4;
% pml = pml + pml.';   pot = pot - 1i*pml;
% clear idx xrange;

[~, ~, dt] = getfft(u0, dx, dy);

qsp = QSParameters({x, y}, dt, pot, zeros(size(u0)));   qsp.q = -1;
% qsp = QSParameters({x, y}, dt, pot);
clear x y r theta;

solver = QuantumSolver2D(qsp);

ut0 = u0;
ut1 = zeros(size(u0));
Et = zeros(tmax,1);

% we'll sample our function at (kx, ky); this gets us a nice distribution:
nk = 20;
k = linspace(0, rmax, nk) / rmax * n/2;
phi = (0:(nk-1))*(pi*(1-sqrt(5)));
kx = round(k.*cos(phi) + n/2);
ky = round(k.*sin(phi) + n/2);

V = zeros(size(ut0));
Ax = zeros(size(ut0));
Ay = zeros(size(ut0));

if want2plot
    fig = mesh(qsp.x, qsp.y, abs(ut0));
    axis([-rmax rmax -rmax rmax 0 .6]);
end

tic;
for t = 1:tmax
    if isempty(qsp.mag)
        ut0 = solver.solve(ut0);
    else
        ut0 = solver.solve(ut0, ut1);
        ut1 = solver.solve(ut1, ut0);
        
        [V, Ax, Ay] = solver.solvePotentials(ut0, ut1);
        solver.updateHamiltonian(V, Ax, Ay);
    end
    
    for i = 1:length(k)
        Et(t) = Et(t) + ut0(kx(i), ky(i));
    end
    
    if want2plot && ~mod(t-1, 1)
        fig.ZData = abs(ut0);
        drawnow;
    end
    
    displayProgress(t-1, tmax);
end

qsr = QSResults(Et, dt);
qsr.spectrum;
drawnow;

clear solver.solve;
clear mex; %#ok<CLMEX>
% finished;

function displayProgress(t, tmax)
    progress_counter = round(0.05 * tmax);
    
    if ~mod(t, progress_counter)
        progress_text = ['Progress: ', num2str(round(t / tmax * 100)), '%'];
        disp([newline, progress_text]);
        toc
    end
end

% function finished
%     for t = 1:5
%         beep; pause(.5);
%     end
% end