function [X,Y,Xu,Yu,Xv,Yv,hx,hy] = build_EulerianGrid(Nx,Ny,La,Lb,Lc,Ld,BC_choice)

if BC_choice==0 
    hx = (Lb-La)/Nx; hy = (Ld-Lc)/Ny;
    x = hx*[1:Nx]-hx+La;
    y = hy*[1:Ny]-hy+Lc;
    xh = x+hx/2;
    yh = y+hy/2;
    
else
    hx = (Lb-La)/(Nx+1); hy = (Ld-Lc)/(Ny+1);
    x = hx*[0:Nx+1] + La;   % including endpoints
    y = hy*[0:Ny+1] + Lc;   % including endpoints
    xh = x(1:end-1)+hx/2;
    yh = y(1:end-1)+hy/2;
end

[X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
X = X';                     % transpose so that X(i,j),Y(i,j) are
Y = Y';                     % coordinates of (i,j) point

if BC_choice == 0
    % Xu = X; Yu = Y;
    % Xv = X; Yv = Y;
    % points where u values are solved
    [Xu,Yu] = meshgrid(x,yh);
    Xu = Xu'; Yu = Yu';
    % points where v values are solved
    [Xv,Yv] = meshgrid(xh,y);
    Xv = Xv'; Yv = Yv';
else
    % points where u values are solved
    if BC_choice==2
        [Xu,Yu] = meshgrid(x(1:end-1),yh);
    else
        [Xu,Yu] = meshgrid(x,yh);
    end
    Xu = Xu'; Yu = Yu';
    % points where v values are solved
    [Xv,Yv] = meshgrid(xh,y);
    Xv = Xv'; Yv = Yv';
end

end