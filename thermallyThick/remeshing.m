% Remeshing
function [nx_new, dx_new] = remeshing(nx_old, delta)

global dx_i

% Design objective:
% - Maintain uniform grid with spacing close to initial value
dx0 = dx_i;   % Initial grid cell size

% Update the number of cells
nx_new = min(round(delta/dx0),nx_old);   % nx_new <= nx_old
nx_new = max(nx_new,5);                  % Keep minimum of 5 grid cells

% Update grid cell size
dx_new = (delta/nx_new)*ones(nx_new,1);

end


 