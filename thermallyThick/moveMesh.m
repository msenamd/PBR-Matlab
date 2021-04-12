function [xRight, xCenter, dx] = moveMesh(nx_old,dV,geometry,A_rectangle,L_cylinder)

%reconstructing the mesh (Note: dx~dr , xRight~r , xCenter~r_c)
if geometry=='rectangle'    
    dx=dV/A_rectangle;
    xRight(1) = dx(1);
    xCenter(1) = 0.5*xRight(1);
    for i = 2:nx_old
    xRight(i) = dx(i) + xRight(i-1);
    xCenter(i) = 0.5*( xRight(i-1) + xRight(i) );
    end
          
elseif geometry=='cylinder'
    dx(1)= (dV(1)/L_cylinder/pi)^(1/2);
    xRight(1) = dx(1);
    xCenter(1) = 0.5*xRight(1);
    for i=2:nx_old
        dx(i)= ((dV(i)/L_cylinder+pi*xRight(i-1)^2)/pi)^0.5-xRight(i-1);
        xRight(i) = dx(i) + xRight(i-1);
        xCenter(i) = 0.5*( xRight(i-1) + xRight(i) );
    end
    
elseif geometry=='sphere'
    dx(1)= (3*dV(1)/4/pi)^(1/3);
    xRight(1) = dx(1);
    xCenter(1) = 0.5*xRight(1);
    for i=2:nx_old
        dx(i)= (3*dV(i)/4/pi+xRight(i-1)^3)^(1/3)-xRight(i-1);
        xRight(i) = dx(i) + xRight(i-1);
        xCenter(i) = 0.5*( xRight(i-1) + xRight(i) );
    end


end
