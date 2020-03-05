classdef (Abstract) QuantumSolvers < handle
% QuantumSolvers   Interface for classes solving non-relativistic problems
%   in quantum mechanics.
%
%   See also: QuantumSolver1D, QuantumSolver2D.
    
    properties (Abstract, SetAccess = private, GetAccess = public)
        qs_parameters
        magnetic_field
    end
    
    methods (Abstract)
        solve(wf, qsp)
    end
    
    methods (Abstract, Access = protected, Static)
        setup(qsp)
    end
    
    methods (Static)
        % Checks if your time step is small enough so that aliasing effects
        % are not present; or at least minimal.
        function [alias, dt] = willalias(u, qsp)
            if isvector(u)
                [~, ~, dt] = CommonFunctions.getfft(u, qsp.dx, qsp.m);
            elseif ismatrix(u)
                [~, ~, dt] = CommonFunctions.getfft(u, qsp.dx, qsp.dy, qsp.m);
            end
            
            alias = dt < qsp.dt;
        end
    end
end