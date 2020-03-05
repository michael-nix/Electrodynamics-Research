% isa(QuantumSolver2D, 'handle') == true
classdef QuantumSolver2D < QuantumSolvers 
    
    properties (SetAccess = private, GetAccess = public)
        qs_parameters
        magnetic_field
        lhs1, lhs2, rhs1, rhs2
    end
    
    methods
        function solver = QuantumSolver2D(qsp)
            solver.qs_parameters = qsp;
            
            [solver.lhs1, solver.rhs1, ...
                solver.lhs2, solver.rhs2] = solver.setup(qsp);
            
            if ~isempty(qsp.mag)
                m = qsp.m;   h_bar = qsp.h_bar;
                
                if isa(qsp.mag, 'function_handle')
                    solver.magnetic_field = 1i * h_bar / 2 / m * qsp.mag(qsp.x);
                    
                elseif isnumeric(qsp.mag)
                    solver.magnetic_field = 1i * h_bar / 2 / m * qsp.mag;
                    
                else
                    error('Magnetic Field improperly defined. Must be matrix or function handle.');
                end
            end
        end
        
        function ut = solve(solver, u0, u1)
            import CommonFunctions.*;
            
            % --- Step One:
            uts = solver.rhs1(:,:,2).*u0;
            uts(2:end,:) = uts(2:end,:) + solver.rhs1(1:end-1,:,1).*u0(1:end-1,:);
            uts(1:end-1,:) = uts(1:end-1,:) + solver.rhs1(2:end,:,3).*u0(2:end,:);
            
            if nargin == 3
                if ~isempty(solver.magnetic_field)
                    uts = uts + solver.magnetic_field .* u1;
                    if isreal(uts)
                        uts = complex(uts); 
                    end
                    
                else
                    error('No magnetic field defined.');
                end
            end
            
%             uts = tdma(solver.lhs1(:,:,2).', solver.lhs1(:,:,1).', solver.lhs1(:,:,3).', uts.').';
            uts = tdma(solver.lhs1(:,:,2), solver.lhs1(:,:,1), solver.lhs1(:,:,3), uts.').';
            
            % --- Step Two:
            ut = solver.rhs2(:,:,2).*uts;
            ut(:,2:end) = ut(:,2:end) + solver.rhs2(:,1:end-1,1).*uts(:,1:end-1);
            ut(:,1:end-1) = ut(:,1:end-1) + solver.rhs2(:,2:end,3).*uts(:,2:end);
            
            if nargin == 3
                if ~isempty(solver.magnetic_field)
                    ut = ut + solver.magnetic_field .* u1;
                    if isreal(ut)
                        ut = complex(ut); 
                    end
                else
                    error('No magnetic field defined.');
                end
            end
            
            ut = tdma(solver.lhs2(:,:,2), solver.lhs2(:,:,1), solver.lhs2(:,:,3), ut);
        end
        
        function [V_next, Ax_next, Ay_next] = solvePotentials(solver, u0, u1)
            dr = solver.qs_parameters.dx;
            c = solver.qs_parameters.c;
            dt = solver.qs_parameters.dt;
                        
            persistent sigmax sigmay s_xtimesy Cvx Bvx Cvy Bvy Cu Bu;
            if isempty(sigmax)
                x = solver.qs_parameters.x;   
                
                [sigmax, sigmay] = solver.setupPML(x, dr);
                s_xtimesy = sparse(sigmax.*sigmay / c^2);
                
                Cvx = 1./(1 + dt/2 * sigmax);   Bvx = 1 - dt/2 * sigmax;
                Cvy = 1./(1 + dt/2 * sigmay);   Bvy = 1 - dt/2 * sigmay;
                Cu = 1./(1 + dt/2 * (sigmax + sigmay));   
                Bu = 1 - dt/2 * (sigmax + sigmay);
            end
            
            % ---
            persistent V_now vx vy psi
            if isempty(V_now)
                m = size(u0);
                V_now = zeros(m);   
                psi = zeros(m);
                vx = zeros(m);   vy = zeros(m);
            end
            
            dudx = solver.gradient_x(V_now, dr);
            dudy = solver.gradient_y(V_now, dr);
            vx = Cvx.*(Bvx.*vx + dt*dudx);
            vy = Cvy.*(Bvy.*vy + dt*dudy);
            
            dudx = solver.gradient_x(vx, dr);
            dudy = solver.gradient_y(vy, dr);
            psi = psi + dt * (sigmay.*dudx + sigmax.*dudy - s_xtimesy.*V_now + u0.*conj(u0) + u1.*conj(u1));
            
            V_now = Cu.*(Bu.*V_now + dt*c^2 * (dudx + dudy + psi));
            V_next = V_now;
            
            % ---
            persistent Ax_now Axvx Axvy Axpsi
            if isempty(Ax_now)
                m = size(u0);
                Ax_now = zeros(m);   
                Axpsi = zeros(m);
                Axvx = zeros(m);   Axvy = zeros(m);
            end
            
            dudx = solver.gradient_x(Ax_now, dr);
            dudy = solver.gradient_y(Ax_now, dr);
            Axvx = Cvx.*(Bvx.*Axvx + dt*dudx);
            Axvy = Cvy.*(Bvy.*Axvy + dt*dudy);
            
            dudx = solver.gradient_x(Axvx, dr);
            dudy = solver.gradient_y(Axvy, dr);
            Axpsi = Axpsi + dt * (sigmay.*dudx + sigmax.*dudy - s_xtimesy.*Ax_now + conj(u0).*u1 + conj(u1).*u0);
            
            Ax_now = Cu.*(Bu.*Ax_now + dt*c^2 * (dudx + dudy + Axpsi));
            Ax_next = Ax_now;

            % ---
            persistent Ay_now Ayvx Ayvy Aypsi
            if isempty(Ay_now)
                m = size(u0);
                Ay_now = zeros(m);   
                Aypsi = zeros(m);
                Ayvx = zeros(m);   Ayvy = zeros(m);
            end
            
            dudx = solver.gradient_x(Ay_now, dr);
            dudy = solver.gradient_y(Ay_now, dr);
            Ayvx = Cvx.*(Bvx.*Ayvx + dt*dudx);
            Ayvy = Cvy.*(Bvy.*Ayvy + dt*dudy);
            
            dudx = solver.gradient_x(Ayvx, dr);
            dudy = solver.gradient_y(Ayvy, dr);
            Aypsi = Aypsi + dt * (sigmay.*dudx + sigmax.*dudy - s_xtimesy.*Ay_now - 1i*conj(u0).*u1 + 1i*conj(u1).*u0);
            
            Ay_now = Cu.*(Bu.*Ay_now + dt*c^2 * (dudx + dudy + Aypsi));
            Ay_next = Ay_now;
            
            % ---
            dudx = solver.gradient_x(Ay_next, dr);
            dudy = solver.gradient_y(Ax_next, dr);
            solver.magnetic_field = dudx - dudy;
        end
        
        function updateHamiltonian(solver, V, Ax, Ay)
            qsp = solver.qs_parameters;
            
            dt = qsp.dt / 2;    dx = qsp.dx;    dy = qsp.dy;    n = qsp.n;
            h_bar = qsp.h_bar;    m = qsp.m;    q = qsp.q;
            nx = n(2);    ny = n(1);
            
            a = 1i / h_bar / 2 / m;
            
            dAxdx = solver.gradient_x(Ax, dx);
            dAydy = solver.gradient_y(Ay, dy);
            d = -1i * q * V / h_bar ...
                + q / m / 2 * ((dAxdx + dAydy) - 1i * q * (Ax.^2 + Ay.^2));
            
            b = q * h_bar * Ax / 2 / m;
            c = q * h_bar * Ay / 2 / m;
            
            % --- Step 1:
            solver.lhs1(:,:,1) = (-a/dx^2 + b/dx/2).';         % diag - 1
            solver.lhs1(:,:,2) = (1/dt + 2*a/dx^2 - d/2).';    % diag
            solver.lhs1(:,:,3) = (-a/dx^2 - b/dx/2).';         % diag + 1
            
            solver.rhs1(:,:,1) = a/dy^2 - c/dy/2;          % diag - 1
            solver.rhs1(:,:,2) = 1/dt - 2*a/dy^2 + d/2;    % diag
            solver.rhs1(:,:,3) = a/dy^2 + c/dy/2;          % diag + 1
            
            % --- Step 2:
            solver.lhs2(:,:,1) = -a/dy^2 + c/dy/2;         % diag - 1
            solver.lhs2(:,:,2) = 1/dt + 2*a/dy^2 - d/2;    % diag
            solver.lhs2(:,:,3) = -a/dy^2 - c/dy/2;         % diag + 1
            
            solver.rhs2 = zeros(ny, nx, 3);
            solver.rhs2(:,:,1) = a/dx^2 - b/dx/2;         % diag - 1
            solver.rhs2(:,:,2) = 1/dt - 2*a/dx^2 + d/2;   % diag
            solver.rhs2(:,:,3) = a/dx^2 + b/dx/2;         % diag + 1
            
        end
    end
    
    methods (Access = protected, Static)
        function [lhs1, rhs1, lhs2, rhs2] = setup(qsp)
            dt = qsp.dt / 2; dx = qsp.dx; dy = qsp.dy; n = qsp.n;
            h_bar = qsp.h_bar; m = qsp.m; q = qsp.q;
            nx = n(2); ny = n(1);
            
            if isa(qsp.pot, 'function_handle')
                V = qsp.pot(sqrt(qsp.x.^2 + qsp.y.^2));
            elseif isnumeric(qsp.pot)
                V = qsp.pot;
            else
                error('Scalar Potential improperly defined.  Must be numeric matrix or function handle.');
            end
            
            % ut = a (uxx + uyy) + b ux + c uy + d u
            a = 1i / h_bar / 2 / m;
            d = -1i * q * V / h_bar;
            
            b = 0;
            c = 0;
            
            % --- Step 1:
            lhs1 = zeros(ny, nx, 3);
%             lhs1(:,:,1) = -a/dx^2 + b/dx/2;         % diag - 1
%             lhs1(:,:,2) = 1/dt + 2*a/dx^2 - d/2;    % diag
%             lhs1(:,:,3) = -a/dx^2 - b/dx/2;         % diag + 1
            lhs1(:,:,1) = (-a/dx^2 + b/dx/2).';         % diag - 1
            lhs1(:,:,2) = (1/dt + 2*a/dx^2 - d/2).';    % diag
            lhs1(:,:,3) = (-a/dx^2 - b/dx/2).';         % diag + 1
            
            rhs1 = zeros(ny, nx, 3);
            rhs1(:,:,1) = a/dy^2 - c/dy/2;          % diag - 1
            rhs1(:,:,2) = 1/dt - 2*a/dy^2 + d/2;    % diag
            rhs1(:,:,3) = a/dy^2 + c/dy/2;          % diag + 1
            
            % --- Step 2:
            lhs2 = zeros(ny, nx, 3);
            lhs2(:,:,1) = -a/dy^2 + c/dy/2;         % diag - 1
            lhs2(:,:,2) = 1/dt + 2*a/dy^2 - d/2;    % diag
            lhs2(:,:,3) = -a/dy^2 - c/dy/2;         % diag + 1
            
            rhs2 = zeros(ny, nx, 3);
            rhs2(:,:,1) = a/dx^2 - b/dx/2;         % diag - 1
            rhs2(:,:,2) = 1/dt - 2*a/dx^2 + d/2;   % diag
            rhs2(:,:,3) = a/dx^2 + b/dx/2;         % diag + 1
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
        
        function g = divergence(fx, fy, dx)
            g = gradient_x(fx, dx) + gradient_y(fy, dx);
        end
        
        function [sigmax, sigmay] = setupPML(x, dx)
            PMLwidth = 15;   m = 4;
            
            PMLwidth = PMLwidth * dx;
            sigmax_max = log(1e4)*(m+1) / PMLwidth^(m+1);
            PMLwidth = max(max(x)) - PMLwidth;
            
            sigmax = zeros(size(x));
            idx = abs(x) > PMLwidth;
            sigmax(idx) = sigmax_max * (abs(x(idx)) - PMLwidth).^m;
            sigmay = sigmax.';
            
            sigmax = sparse(sigmax);
            sigmay = sparse(sigmay);
        end
    end
end