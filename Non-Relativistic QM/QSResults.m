classdef QSResults < handle
% QSRESULTS  Object that parses results and provides easy plots of results
% from simulations run with QuantumSolvers using QSParameters objects.
%   
%   qsr = QSResults(uts, dt) takes in the time signal, uts, generated by a
%   QuantumSolver, and its time step, dt, to calculate the energy spectrum
%   and mode energies present.
%   
%   QSResults methods:
%       spectrum    - plots the energy spectrum.
%       modes       - plots the modes and their energies.
%       timesignal  - plots the time signal used to create this object.
%   
%   See also QSParameters, QuantumSolver1D, MainQuantumSolverUI.
    properties (GetAccess = public, SetAccess = private)
        uts = 0; dt = 0;
        energy = 0; amplitude = 0;
        mode_energies = 0;
        qs_parameters = 0;
    end
    
    properties (GetAccess = private, SetAccess = private)
        peaks
    end
    
    methods
        function qsr = QSResults(uts, dt, qsp)
            if nargin == 0
                return;
            end
            
            if nargin == 3
                qsr.qs_parameters = qsp;
            end
            
            import CommonFunctions.getfft;
            qsr.uts = uts;
            qsr.dt = dt;
            
            [amplitude, energy] = getfft(uts, dt);
            qsr.amplitude = amplitude;
            qsr.energy = -energy;
            
            qsr.peaks = islocalmax(abs(amplitude)) & abs(amplitude) > (0.01*max(abs(amplitude)));
            qsr.mode_energies = flip(qsr.energy(qsr.peaks == 1));
        end
        
        function spectrum(qsr, axes_obj)
            if nargin < 2
                figure;
            end
            
            xmin = qsr.mode_energies(1) - 1;
            xmax = qsr.mode_energies(end) + 1;
            
            hold on;
            if nargin < 2
                plot(qsr.energy, abs(qsr.amplitude));
                plot(qsr.energy(qsr.peaks), abs(qsr.amplitude(qsr.peaks)), 'rx');
            else
                plot(axes_obj, qsr.energy, abs(qsr.amplitude));
                plot(axes_obj, qsr.energy(qsr.peaks), abs(qsr.amplitude(qsr.peaks)), 'rx');
            end
            axis([xmin, xmax, 0, max(abs(qsr.amplitude))*1.1]);
            xlabel('Energy'); ylabel('Amplitude');
            title('Wave Function Energy Spectrum');
        end
        
        function modes(qsr, axes_obj)
            if nargin < 2
                figure;
            end
            
            modes = 1:length(qsr.mode_energies);
            
            if nargin < 2
                plot(modes, qsr.mode_energies, '-o', 'LineWidth', 1.5);
            else
                plot(axes_obj, modes, qsr.mode_energies, '-o', 'LineWidth', 1.5);
            end
            xlabel('Mode Number'); ylabel('Energy'); grid on;
            title('Wave Function Mode Energies');
        end
        
        function timesignal(qsr, axes_obj)
            if nargin < 2
                figure;
            end
            
            t = 0:qsr.dt:(qsr.dt*(length(qsr.uts)-1));
            
            if nargin < 2
                plot(t, abs(qsr.uts));
            else
                plot(axes_obj, t, abs(qsr.uts));
            end
            xlabel('Time'); ylabel('Sum of Amplitudes'); grid on;
            title('Wave Function Measured Over Time');
        end
        
        function qsr = copy(qsr, qsr2)
            qsr.uts = qsr2.uts;
            qsr.dt = qsr2.dt;
            qsr.energy = qsr2.energy;
            qsr.amplitude = qsr2.amplitude;
            qsr.peaks = qsr2.peaks;
            qsr.mode_energies = qsr2.mode_energies;
            qsr.qs_parameters = qsr2.qs_parameters;
        end
    end
    
    methods (Static)
        function qsr = empty
            qsr = QSResults;
        end
    end
end