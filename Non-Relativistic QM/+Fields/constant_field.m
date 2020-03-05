function B = constant_field(x, Bo)
%   CONSTANT_FIELD  Returns a constant magnetic field.
%
%   B = CONSTANT_FIELD(x) returns a constant magnetic field at any or all
%       points given by x.
%   
%   B = CONSTANT_FIELD(x, Bo) returns the value of Bo at any or all points
%       given by x.
%   
%   This function exists to set a standard format for all other fields,
%   mostly applicable for two-dimensional problems with vector potentials.
%
%   Since the use of all fields is as function handles, another way to
%   modify the magnitude of the field would be through
%       mag = @(x) Bo * constant_field(x);
%
%   This will work as long as constant_field.m is in the current path or
%   has been imported with the Fields package.  This applies equally well
%   to all of the potentials.

if nargin == 2
    B = Bo * ones(size(x));
else
    B = 0.4 * ones(size(x));
end