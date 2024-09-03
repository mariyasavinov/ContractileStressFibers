function [u,v] = fluid_2D_solver(XYuv,A,f0,fmyosin,Ffiber_spread,BC_choice)

X = XYuv{1}; %Y = XYuv{2};
% Xu = XYuv{3}; Yu = XYuv{4};
% Xv = XYuv{5}; Yv = XYuv{6};

% Get number of X pts and number of Y pts
if BC_choice==0
    Nx = size(X,1);
    Ny = size(X,2);
else
   Nx = size(X,1)-2;
   Ny = size(X,2)-2;
end

% total force
ftotal = - fmyosin - f0*Ffiber_spread;

% solve for the velocity field
w = A\ftotal;

if BC_choice==0
    u = reshape(w(1:Nx*Ny),Nx,Ny);
    v = reshape(w(Nx*Ny+1:end),Nx,Ny);
elseif BC_choice==1
    u = reshape(w(1:(Nx+2)*(Ny+1)),Nx+2,Ny+1);
    v = reshape(w((Nx+2)*(Ny+1)+1:end),Nx+1,Ny+2);
elseif BC_choice==2
    u = reshape(w(1:(Nx+1)*(Ny+1)),Nx+1,Ny+1);
    v = reshape(w((Nx+1)*(Ny+1)+1:end),Nx+1,Ny+2);
end

end