classdef QSParameters < handle
% QSPARAMETERS  Object that defines the scope of problems for QuantumSolver
%   objects such as QuantumSolver1D or QuantumSolver2D.
%
%   QSParameters objects are defined by their usability by QuantumSolver
%   objects.  As long as a valid range of x-values and a time step,
%   dt, are defined, QSParameters objects can immediately be used to
%   create QuantumSolver objects to simulate simple wavefunctions.
%
%   ps = QSParameters(x, dt) creates a new QSParameters object with
%       default parameters and potential.
%
%   ps = QSParameters(x, dt, pot) creates a new QSParameters object
%       with default parameters and classical electric potential
%       defined by a function handle given by pot.
%
%   ps = QSParameters(x, dt, pot, mag) creates a new QSParameters
%       object with electric potential and one dimensional magnetic
%       field defined by a function handle given by mag.
%
%   Examples
%       dt = 0.1;
%       x = -5:0.1:5;
%       ps = QSParameters(x, dt);
%
%       pot = @harmonic_potential;
%       ps = QSParameters(x, dt, pot);
%
%       magfield = @constant_field;
%       ps = QSParameters(x, dt, function_handle.empty, magfield);
%
%   See also QSResults, QuantumSolver1D, QuantumSolvers.
    properties
        x, dx, y, dy, z, dz, n, dt, Ax
        h_bar = 1; m = 1; c = 1; q = 1;
        pot = @Potentials.zero_potential;
        mag
        wavefunction
        wavefunction_time
        wavefunction_spin
        samples
        tmax = 16384;
    end
    
    methods
        function ps = QSParameters(x, dt, pot, mag, m)
            if nargin < 2
                error('QSParameters:InvalidInput','Not enough input arguments.');
            end
            
            if ~iscell(x)
                if ~isvector(x)
                    error('QSParameters:InvalidInput','For 2D problems x must be a two-element cell array');
                end
                if isrow(x)
                    x = x.';
                end
                
                ps.x = x;
                ps.dx = abs(x(2) - x(1));
                ps.n = length(x);
                ps.dt = dt;
                ps.samples = floor(ps.n ./ ((1:5) + 1));
            else
                if length(x) ~= 2
                    error('QSParameters:InvalidInput','For 2D problems x must be a two-element cell array');
                end
                ps.x = x{1};
                ps.y = x{2};
                
                ps.dx = abs(ps.x(1,2) - ps.x(1));
                ps.dy = abs(ps.y(2) - ps.y(1));
                ps.n = size(ps.x);
                ps.dt = dt;
            end
            
            if nargin > 2
                if isa(pot, 'function_handle')
                    if ~isempty(pot)
                        ps.pot = pot;
                    else
                        ps.pot = @zero_potential;
                    end
                elseif isnumeric(pot)
                    ps.pot = pot;
                else
                    error('QSParameters:InvalidInput','Invalid potential: must be numeric array or non-empty function handle.');
                end
            end
            
            if nargin > 3
                if isa(mag, 'function_handle') && ~isempty(mag)
                    ps.mag = mag;
                elseif isnumeric(mag)
                    ps.mag = mag;
                else
                    error('QSParameters:InvalidInput','Invalid field: must be numeric array or non-empty function handle.');
                end
            end
            
            if nargin > 4
                ps.m = m;
            end
        end
    end
end