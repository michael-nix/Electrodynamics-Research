function [U, p, dt] = getfft(u, dx, varargin)
% GETFFT  Helper function for Fourier Transforms and to check for aliasing.
%   It's like FFTSHIFT but also formats the momentum properly.  Also, it
%   makes sure everything is a power of two and returns things as such.
%   
%   Available function signatures:
%   
%   U = getfft(u) 1D or 2D
%   
%   [U, p] = getfft(u, dx) 1D only
%   
%   [U, p, dt] = getfft(u, dx, m) 1D only
%   
%   [U, p] = getfft(u, dx, dy) 2D only
%   
%   [U, p, dt] = getfft(u, dx, dy, m) 2D only
    
    InfOrNaN = find(~isfinite(u), 1);
    if ~isempty(InfOrNaN)
        error(['Inf or NaN will break your results.  You have at least one' ...
            ', at index ', num2str(InfOrNaN), ' of your input function.']);
    end
    
    if ~ismatrix(u)
        error('Only supports 1D & 2D input.');
    end
    
    if nargin == 3
        if isvector(u)
            m = varargin{1};
        else
            dy = varargin{1};
        end
    elseif nargin == 4
        dy = varargin{1};
        m = varargin{2};
    end
    
    if isvector(u)
        n = 2^nextpow2(length(u));
        U = fft(u, n)/n;
        
        U = fftshift(U);
        
        if nargout > 1
            p = ((0:n-1) - n/2) / n * (2*pi/dx);
            if iscolumn(u)
                p = p.';
            end
            
            if nargout > 2
                p_abs = abs(p);
            end
        end
    elseif ismatrix(u)
        n = 2.^nextpow2(size(u));
        U = fft2(u, n(1), n(2))/prod(n);
        U = fftshift(U);
        
        if nargout > 1
            if ~exist('dx', 'var')
                error('Not enough input arguments. Variable ''dx'' undefined.');
            end
            if ~exist('dy', 'var')
                error('Not enough input arguments. Variable ''dy'' undefined.');
            end
            
            px = ((0:n(1)-1) - n(1)/2) / n(1) * (2*pi/dx);
            py = ((0:n(2)-1) - n(2)/2) / n(2) * (2*pi/dy);
            [p.x, p.y] = meshgrid(px, py);
            
            if nargout > 2
                p_abs = sqrt(p.x.^2 + p.y.^2);
            end
        end
    end
    
    if nargout > 2
        if ~exist('m', 'var')
            % warning('Not enough input arguments. Default value of m = 1.');
            m = 1;
        end
        
        V = abs(U);
        V = V / max(V(:));
        
        pmax = 2 * max(p_abs(V > exp(-5)));
        if isempty(pmax)
            pmax = max(p_abs(:));
        end
        
        emax = pmax^2 / 2 / m;
        dt = 1/emax;
    end