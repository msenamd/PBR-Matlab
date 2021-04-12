% Interpolate solution on new mesh
function [temp, x_m, x_vs, x_c] = interpolateOnNewMesh(nx_new,xCenter, ...
                     nx_old,xCenter_old,temp_old,x_m_old,x_vs_old,x_c_old)
for i=1:nx_new
    
    ii = 0;
    for j=1:(nx_old-1)
        if( ( (xCenter_old(j)-xCenter(i))   <= 0 ) ...
          & ( (xCenter(i)-xCenter_old(j+1)) < 0  ) )
            ii = j;
        end
    end
    if( (ii == 0) & (xCenter(i) < xCenter_old(1)) )
        ii = 1;
    end
    
    if( (ii == 0) & (xCenter_old(nx_old) <= xCenter(i)) )
        ii = nx_old-1;
    end
    
    if(ii > 0) 
        temp(i) = temp_old(ii) ...
     +(temp_old(ii+1)-temp_old(ii))*(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));
        x_m(i) = x_m_old(ii) ...
       +(x_m_old(ii+1)-x_m_old(ii))*(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));
        x_vs(i) = x_vs_old(ii) ...
     +(x_vs_old(ii+1)-x_vs_old(ii))*(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));
        x_c(i) = x_c_old(ii) ...
       +(x_c_old(ii+1)-x_c_old(ii))*(xCenter(i)       -xCenter_old(ii)) ...
                                   /(xCenter_old(ii+1)-xCenter_old(ii));
    end
    
    if(ii == 0)
        fprintf(' *** Error in interpolation routine, i = %g \n',i);
        return;
    end

end


 