% Computational grid
%
% Definitions
% - nx:      number of cells
% - dx:      array containing sizes of cells
% - xRight:  array containing coordinates of right edges of cells
% - xCenter: array containing coordinates of cell centers
% - xLeft1:  coordinate of left edge of computational domain
%
% We note delta the coordinate of right edge of computational domain
% - delta may be time-dependent (for cases with volume change)
% - By construction: delta = sum_{i=1}^{i=nx} dx(i)
%
% |                                             (exposed surface)
% |                <-- dx(i) -->                x = delta = xRight(nx)
% |----------| ... |-----o-----| ... |----------| 
% |                      |     |
% |                      o x = xCenter(i)
% |                            | x = xRight(i)
% x = xLeft1 = 0
% (back surface)
%

%     (Note: dx~dr , xRight~r , xCenter~r_c)

function [xRight, xCenter, xLeft1, dV] = mesh(nx,dx,geometry...
                                              ,A_rectangle,L_cylinder)

xRight(1) = dx(1);
for i = 2:nx
    xRight(i) = dx(i) + xRight(i-1);
end

xLeft1 = 0; 
xCenter(1) = 0.5*( xLeft1 + xRight(1) );
for i = 2:nx
    xCenter(i) = 0.5*( xRight(i-1) + xRight(i) );
end


if geometry=="rectangle"    
   dV=dx*A_rectangle;
elseif geometry=="cylinder"
    dV(1)=pi*dx(1)*(xRight(1))*L_cylinder;
    for i=2:nx
    dV(i)=pi*dx(i)*(xRight(i)+xRight(i-1))*L_cylinder;
    end
elseif geometry=="sphere"
    dV(1)=4/3*pi*dx(1)*(xRight(1)^2);
    for i=2:nx
    dV(i)=4/3*pi*dx(i)*...
        (xRight(i)^2+xRight(i)*xRight(i-1)+xRight(i-1)^2);
    end   
else
    fprintf('Incorrect geometry type \n');
    fprintf('Available types: rectangle / sphere / cylinder \n');
    return;
end



end




 