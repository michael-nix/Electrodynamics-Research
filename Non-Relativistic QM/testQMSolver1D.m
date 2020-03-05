clear; import CommonFunctions.*; import Potentials.*; import Fields.*; 
% TESTQMSOLVER1D  Script to test QuantumSolver1D and run some experiments.
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
%   See also QSPARAMETERS, QUANTUMSOLVER2D, QUANTUMSOLVERS.

pot = @electric_potential; mag = @constant_field;

n = 1024; tmax = 16384;

if isequal(pot, @zero_potential)
    xmin = -2.5; xmax = 2.5;
    x0 = 0; v0 = 2;
    
elseif isequal(pot, @harmonic_potential)
    xmin = -10; xmax = 10;
    x0 = 0; v0 = 2;
    
elseif isequal(pot, @gaussian_potential)
    xmin = -10; xmax = 10;
    x0 = 0; v0 = 0;
    
elseif isequal(pot, @harmonic_gaussian_potential)
    xmin = -10; xmax = 10;
    x0 = 2; v0 = 2;
    
elseif isequal(pot, @electric_potential)
    xmin = -200; xmax = 200; n = 4096;
    x0 = 5; v0 = 2; tmax = 2^20;
    
elseif isequal(pot, @morse_potential)
    xmin = -15; xmax = 15; De = 10;
    x0 = 0; v0 = 1;
    
else
    xmin = -5; xmax = 5;
    x0 = 0; v0 = 0;
end

n = 2^nextpow2(n); 
x = linspace(xmin, xmax, n).';

while ~isempty(find(x == 0, 1))
    n = 2^(nextpow2(n) + 1);
    if n > 2^12
        error('x has a zero and n is too big.');
    end
    x = linspace(xmin, xmax, n).';
end
dx = x(2) - x(1);

u0 = initeven(x, x0, -v0, 1/sqrt(2)) + initodd(x, -x0, v0, 1/sqrt(2));
u0 = u0 / norm(u0 * sqrt(dx));
[U0, p, dt] = getfft(u0, dx);

qsp = QSParameters(x, dt, pot);
% qsp = QSParameters(x, dt, pot, mag);

solver = QuantumSolver1D(qsp);

ut0 = u0;
Et = zeros(1, tmax);

if isempty(qsp.mag)
    [lhs, rhs] = solver.getSparseInput;
else
    [lhs, rhs, mag] = solver.getSparseInput;
    ut1 = zeros(n, 1);
end

if isempty(qsp.samples)
    samples = round(n./((1:5)+1));
else
    samples = qsp.samples;
end

for t = 1:tmax
    Et(1, t) = sum(ut0(samples));
    
    if isempty(qsp.mag)                     % Alternate methods:
        ut0 = lhs\(rhs*ut0);                % ut0 = solver.solve(ut0);
    else
        ut0 = lhs\(rhs*ut0 + mag.*ut1);     % ut0 = solver.solve(ut0, ut1);
        ut1 = lhs\(rhs*ut1 + mag.*ut0);     % ut1 = solver.solve(ut1, ut0);
    end
    
%     if ~mod(t, 10)
%         plot(x, abs(ut0));
%         drawnow;
%     end
end

if isempty(qsp.mag)
    [Ef, E] = getfft(Et, dt);
    
    dE = E(2) - E(1);
    peaklocs = peakfinder(abs(Ef).^2, 1/tmax/3);
    peaklocs = peaklocs / tmax * (E(end)-E(1)-(dE)) - E(end) - dE;
else
    [Ef, E] = getfft(Et, dt);
    plot(E, abs(Ef));
    axis([-5, 0, 0, 0.5]);
end

if isempty(qsp.mag)
    if isequal(qsp.pot, @zero_potential)
        energy = -peaklocs(end:-1:1);
        elvl = 1:length(energy);
        
        figure;plot(elvl,energy,'ko',elvl,elvl.^2*pi^2/2/(xmax-xmin)^2,'kx');
        axis([elvl(1)-1, elvl(end)+1, 0, max(abs(energy))+1]);
        xlabel('Energy Level (n)'); ylabel('Energy (hartree)');
        
        error = abs(energy-elvl.^2*pi^2/2/(xmax-xmin)^2)./energy*100;
        figure;plot(elvl,error,'kx');axis([0.9 max(elvl)+0.1 0 max(error)*1.05]);
        xlabel('Energy Level (n)'); ylabel('|Error| (%)');
        
    elseif isequal(qsp.pot, @harmonic_potential)
        energy = -peaklocs(end:-1:1);
        elvl = 0:length(energy)-1;
        
        figure; plot(elvl, energy, 'ko', elvl, (2*(elvl)+1)/2, 'kx');
        axis([elvl(1)-1, elvl(end)+1, 0, max(abs(energy))+1]);
        xlabel('Energy Level (n)'); ylabel('Energy (hartree)');
        
        error = abs(energy-(2*(elvl)+1)/2)./energy*100;
        figure; plot(elvl, error, 'kx'); axis tight;
        xlabel('Energy Level (n)'); ylabel('|Error| (%)');
        
    elseif isequal(qsp.pot, @gaussian_potential)
        energy = peaklocs(peaklocs > 0);
        energy = energy(end:-1:1);
        enum = [5.41, 4.28, 3.25, 2.3, 1.22, 0.238];
        
        m = min([length(energy), length(enum)]);
        energy = energy(1:m);
        enum = enum(1:m);
        elvl = 1:m;
        
        figure;plot(elvl,energy,'ko',elvl,enum,'kx');
        axis([0.9 max(elvl)+0.1 0 6]);
        xlabel('Energy Level (n)');ylabel('Energy (hartree)');
        
        error = abs(energy-enum)./energy*100;
        figure;plot(elvl,error,'kx'); axis tight;
        xlabel('Energy Level (n)'); ylabel('|Error| (%)');
        
    elseif isequal(qsp.pot, @morse_potential)
        energy = -peaklocs(end:-1:1);
        elvl = 0:length(energy)-1;
        a = sqrt(1/2/De); nu0 = a*sqrt(2*De);
        etheory = nu0*(elvl+1/2)-(nu0*(elvl+1/2)).^2/4/De;
        
        figure;plot(elvl,energy,'ko',elvl,etheory,'kx');
        axis([elvl(1)-0.1, elvl(end)+0.1, 0, abs(max(energy))+0.1]);
        xlabel('Energy Level (n)');ylabel('Energy (hartree)');
        
        error = abs(energy-etheory)./energy*100;
        figure;plot(elvl,error,'kx'); axis tight;
        xlabel('Energy Level (n)');ylabel('|Error| (%)');
    end
end

qsr = QSResults(Et, dt);