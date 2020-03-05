classdef QuantumSolver1D < QuantumSolvers
% QuantumSolver1D  Solver for one dimensional quantum systems.
%   
%   Crank-Nicolson method for solving Schrodinger's and Pauli's equations
%   in 1D and taking advantage of the Thomas algorithm for solving
%   tri-diagonal matrices.  This method is unconditionally stable
%   regardless of time step, though aliasing effects may occur if the time
%   step and spatial resolution are not properly determined.
% 
%   solver = QuantumSolver1D(QSParameters) creates a new solver object with
%       parameters determined by the given QSParameters object.
%
%   ut_1 = solver.solve(ut_0) finds the QuantumSolver ut_1 at the next time
%       step given the input QuantumSolver ut_0 and the parameters 
%       determined by QSParameters object used to create the solver.
%
%   ut_up_1 = solver.solve(ut_up_0, ut_down_0) finds the wave function at
%       next time step taking into account spin if the magnetic field is
%       defined by the parameter QSParameters.Bx, and given wave functions
%       for both spin up and spin down at the current time.  If no vector
%       potential is defined, this uses a modified version of Pauli's
%       equation.
%
%   Pauli's equation in one dimension can be simulated if the vector
%   potential is also given through the parameter QSParameters.Ax.  Proper
%   definition of the magnetic field given a vector potential is not
%   handled by this object or the QSParameters object, and both vector
%   potential and magnetic field are defined separately in the QSParameters
%   object.
%
%   Example 1
%       dx = 0.1; dt = 0.1;
%       x = -10:dx:10;
%       ut0 = initeven(x);
%       ps = QSParameters(x, dt, @harmonic_potential);
%       solver = QuantumSolver1D(qsp);
%       ut1 = solver.solve(ut0);
%
%   Example 2
%       solver = QuantumSolver1D(qsp);
%       solver = solver.updateMagneticField(@constant_field);
%       ut0 = initeven(x);
%       ut1 = zeros(1,length(x));
%       for t = 1:1024
%           ut0 = solver.solve(ut0, ut1);
%           ut1 = solver.solve(ut1, ut0);
%           plot(x, abs(ut0).^2, x, abs(ut1).^2);
%           title('QuantumSolvers w/ Spin Up & Spin Down');
%           drawnow;
%       end
%   
%   Example 3
%       solver = QuantumSolver1D(qsp);
%       ut0 = initeven(x);
%       [lhs, rhs] = solver.getSparseInput;
%       for t = 1:1024
%           ut0 = lhs\(rhs*ut0);
%       end
%   
%   QuantumSolver1D methods:
%   
%       solve -- solves the Crank-Nicolson problem for a single time step.
%       
%       updateSolver -- re-creates the solver with a new QSParameters
%       object.
%       
%       updateMagneticField -- re-creates the magnetic_field property.
%       
%       getSpareInput -- returns the Crank-Nicolson matrices so that they
%       can be used manually; > 30% speedup!
%       
%       updatePotential -- re-creates the scalar / vector potential and
%       re-creates the necessary Crank-Nicolson matrices.
%   
%   See also QSParameters, testQMSolver1D, QuantumSolvers.
    
    properties (SetAccess = private, GetAccess = public)
        qs_parameters
        cn_matrices
        magnetic_field
    end
    
    methods
        function solver = QuantumSolver1D(qsp)
            solver.qs_parameters = qsp;
            solver.cn_matrices = solver.setup(qsp);
            
            if ~isempty(qsp.mag)
                dt = qsp.dt; dx = qsp.dx; m = qsp.m;
                if isa(qsp.mag, 'function_handle')
                    solver.magnetic_field = 1i * dt * 2 * dx^2 / m * qsp.mag(qsp.x);
                elseif isnumeric(qsp.mag)
                    solver.magnetic_field = 1i * dt * 2 * dx^2 / m * qsp.mag;
                else
                    error('Magnetic Field improperly defined. Must be matrix or function handle.');
                end
            end
            
            lastwarn('');
        end
        
        function ut = solve(solver, u0, u1)
            n = solver.qs_parameters.n;
            
            lhs = solver.cn_matrices(:,1:n);
            rhs = solver.cn_matrices(:,n+1:end) * u0;
            
            if nargin == 3
                if ~isempty(solver.magnetic_field)
                    rhs = rhs + solver.magnetic_field .* u1;
                else
                    error('No magnetic field defined.');
                end
            end
            
            % --- ut(:,1) = lhs(:,:)^(-1) * rhs(:,1)
            ut = lhs \ rhs;
            
            % only really happens with unstable potentials.
            [~, msgID] = lastwarn;
            if ~isempty(msgID) && isequal(msgID,'MATLAB:nearlySingularMatrix')
                error(['MATLAB''s mldivide threw up a warning, ', msgID,...
                    '. Your time step is probably not small enough,' ...
                    ' or there is some other typo in your code, and'...
                    ' you''re about to get crazy results!']);
            end
        end
        
        function solver = updateSolver(solver, qsp)
            solver.cn_matrices = solver.setup(qsp);
        end
        
        function solver = updateMagneticField(solver, mag)
            qsp = solver.parameters;
            dt = qsp.dt; dx = qsp.dx; m = qsp.m;
            
            if isa(mag, 'function_handle')
                qsp.mag = mag;
                solver.parameters = qsp;
                solver.magnetic_field = 1i * dt * 2 * dx^2 / m * mag(qsp.x);
            elseif isnumeric(mag)
                solver.magnetic_field = 1i * dt * 2 * dx^2 / m * mag;
            end
        end
        
        function [lhs, rhs, magfield] = getSparseInput(solver)
            n = solver.qs_parameters.n;
            cn = solver.cn_matrices;
            
            lhs = cn(:, 1:n);
            rhs = cn(:, n+1:end);
            
            if nargout > 2
                magfield = solver.magnetic_field;
                if isempty(magfield)
                    warning('No magnetic field defined, return value empty.');
                end
            end
        end
        
        function solver = updatePotential(solver, V, Ax)
            qsp = solver.qs_parameters;
            n = qsp.n; dx = qsp.dx; dt = qsp.dt; m = qsp.m;
            
            a = 1i * qsp.h_bar / 2 / m;
            c = -1i * V / qsp.h_bar * 2 * dx^2 * dt;
            
            if nargin < 3
                b = zeros(n, 1);
            else
                b = Ax / 2 / m * dt * dx;
            end
            
            A = (-2*dt*a - b);
            B = (4*dx^2 + 4*dt*a - c);
            C = (-2*dt*a + b);
            lhs = spdiags([C, B, A], -1:1, n, n);

            A = (2*dt*a + b);
            B = (4*dx^2 - 4*dt*a + c);
            C = (2*dt*a - b);
            rhs = spdiags([C, B, A], -1:1, n, n);
            
            solver.cn_matrices = [lhs, rhs];
        end
    end
    
    methods (Access = protected, Static)
        function input = setup(qsp)
            dt = qsp.dt; dx = qsp.dx; n = qsp.n;
            h_bar = qsp.h_bar; m = qsp.m;
            
            if isa(qsp.pot, 'function_handle')
                V = qsp.pot(qsp.x);
            elseif isnumeric(qsp.pot)
                V = qsp.pot;
            else
                error('Scalar Potential improperly defined.  Must be matrix or function handle.');
            end
            
            a = 1i * h_bar / 2 / m;
            c = -1i * V / h_bar * 2 * dx^2 * dt;
            
            if isempty(qsp.Ax)
                b = zeros(n, 1);
            elseif isa(qsp.Ax, 'function_handle')
                b = qsp.Ax(qsp.x) / 2 / m * dt * dx;
            elseif isnumeric(ps.Ax)
                b = qsp.Ax / 2 / m * dt * dx;
            end
            
            A = (-2*dt*a - b);                % u^{t+1}_{n+1}
            B = (4*dx^2 + 4*dt*a - c);        % u^{t+1}_n
            C = (-2*dt*a + b);                % u^{t+1}_{n-1}
            lhs = spdiags([C, B, A], -1:1, n, n);

            A = (2*dt*a + b);
            B = (4*dx^2 - 4*dt*a + c);
            C = (2*dt*a - b);
            rhs = spdiags([C, B, A], -1:1, n, n);
            
            issingular = isinf(condest(rhs)) || isinf(condest(lhs));
            if issingular
                error('QuantumSolver1D:setup:InvalidInput', ...
                    ['Something somewhere has gone to infinity and the problem', ...
                    ' we''ll need to solve has become singular or nearly singular.', ...
                    '  To fix this, try changing (usually by increasing) the spatial', ...
                    ' resolution, i.e. Number of Grid Point, n.']);
            end
            
            input = [lhs, rhs];
        end
    end
end