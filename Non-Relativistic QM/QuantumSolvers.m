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
        
        function g = gradient_x(f, dx)
            g = zeros(size(f));
            
            g(:,1) = (f(:,2) - f(:,1)) / dx;
            g(:,end) = (f(:,end) - f(:,end-1)) / dx;
            
            g(:,2:end-1) = (f(:,3:end) - f(:,1:end-2)) / (2*dx);
        end
        
        function g = gradient_y(f, dy)
            g = zeros(size(f));
            
            g(1,:) = (f(2,:) - f(1,:)) / dy;
            g(end,:) = (f(end,:) - f(end-1,:)) / dy;
            
            g(2:end-1,:) = (f(3:end,:) - f(1:end-2,:)) / (2*dy);
        end
        
        function g = divergence_xy(fx, fy, dx, dy)
            g = gradient_x(fx, dx) + gradient_y(fy, dy);
        end
        
        % assumes your x-dimension goes along the columns, i.e. consists of
        % a single row vector, or for 2D, stacked identical row vectors.
        function [sigmax, sigmay] = setupPML(x, dx)
            
            xiscol = iscolumn(x);
            if xiscol
                x = x.';
            end
            
            PMLwidth = 15;   m = 4;
            sigmax_max = log(1e4)*(m+1) / (PMLwidth*dx)^(m+1);
            
            sigmax = sigmax_max*((1:PMLwidth)*dx).^m;
            sigmax = [sigmax(end:-1:1), ...
                      zeros(1, size(x, 2) - 2 * PMLwidth), sigmax];
            
            sigmax = sparse(sigmax);
            sigmax = repmat(sigmax, size(x, 1), 1);
            
            if xiscol
                sigmax = sigmax.';
            end
            
            if nargout > 1
                sigmay = sigmax.';
            end
        end
    end
end
