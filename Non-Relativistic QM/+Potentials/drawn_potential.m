function V = drawn_potential(x, line_obj)

if nargin < 2 || isequal(line_obj.Visible, 'off')
    x1 = [ -1,   0,   1,  2,   3];
    y1 = [-16, -10, -12, -2, -10];
    y_init = 0;
else
    if ~isa(line_obj, 'matlab.graphics.primitive.Line')
        error('CustomPotential:InvalidInput', 'Custom potentials pull data from Line objects drawn on the Setup axes.');
    end
    x1 = line_obj.XData;
    y1 = line_obj.YData;
    
    [x1, idx] = sort(x1);
    y1 = y1(idx);
    
    y_init = line_obj.UserData;
    if isempty(y_init)
        y_init = 0;
    end
end

V(x < -1) = y_init;
for i = 1:length(x1)
    V(x >= x1(i)) = y1(i);
end

if isrow(V)
    V = V.';
end