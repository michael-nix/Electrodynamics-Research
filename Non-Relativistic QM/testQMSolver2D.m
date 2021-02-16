clear all;   import CommonFunctions.* Fields.* Potentials.*; %#ok<CLALL>
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

pot = @harmonic_potential;

want2plot = false;

if ~exist('+CommonFunctions/tdma.mexw64', 'file')
    error('You''ll need to compile the necessary MEX file from tdma.c');
end

% 2^14 @ 128 = < 40 s.
% 2^14 @ 256 = < 270 s.
% 2^14 @ 512 = < 1100 s.
% 2^14 @ 1024 = < 4450 s.
tmax = 2^14;    n = 256;

rmax = 10;
x = linspace(-rmax, rmax, n);
y = linspace(-rmax, rmax, n);
dx = abs(x(2) - x(1));
dy = abs(y(2) - y(1));
[x, y] = meshgrid(x, y);

% r = sqrt(x.^2+y.^2);
% theta = atan2(y, x);
sig = 1;

u0 = initeven(x, 0, 0, 2).*initeven(y, 0, -2, 2);
u0 = u0 / sqrt(sum(abs(u0).^2 * dx * dy, 'all'));

[~, ~, dt] = getfft(u0, dx, dy);

% qsp = QSParameters({x, y}, dt, pot, zeros(size(u0)));   qsp.q = -100;
qsp = QSParameters({x, y}, dt, pot);
clear x y r theta;

solver = QuantumSolver2D(qsp);

ut0 = u0;
ut1 = 0.1*u0; % zeros(size(u0));
Et = zeros(tmax, 1);

% we'll sample our function at (kx, ky); this gets us a nice distribution:
nk = 20;
k = linspace(0, rmax, nk) / rmax * n/2;
phi = (0:(nk-1)) * (pi*(1 - sqrt(5)));
kx = round(k.*cos(phi) + n/2);
ky = round(k.*sin(phi) + n/2);

if ~isempty(qsp.mag)
    V = zeros(size(ut0));
    Ax = zeros(size(ut0));
    Ay = zeros(size(ut0));
end

if want2plot
    fig = mesh(qsp.x, qsp.y, abs(ut0)); %#ok
    axis([-rmax rmax -rmax rmax 0 0.6]);
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
    
    if want2plot
        fig.ZData = abs(ut0); %#ok
        drawnow;
    end
    
    displayProgress(t-1, tmax);
end

qsr = QSResults(Et, dt);
qsr.spectrum;
drawnow;

clear solver.solve mex; %#ok<CLMEX>

function displayProgress(t, tmax)
    progress_counter = round(0.05 * tmax);
    
    if ~mod(t, progress_counter)
        progress_text = ['Progress: ', num2str(round(t / tmax * 100)), '%'];
        disp([newline, progress_text]);
        toc
    end
end
