function f = myosinForce(XYuv,k,m_func,dmdx_func,dmdy_func,BC_choice)

X = XYuv{1}; Y = XYuv{2};
Xu = XYuv{3}; Yu = XYuv{4};
Xv = XYuv{5}; Yv = XYuv{6};

% Get number of X pts and number of Y pts
if BC_choice==0
    Nx = size(X,1);
    Ny = size(X,2);
else
   Nx = size(X,1)-2;
   Ny = size(X,2)-2;
end


if BC_choice==0
    f1 = k*dmdx_func(Xu,Yu);
    f2 = k*dmdy_func(Xv,Yv);
    f = [reshape(f1,Nx*Ny,1); reshape(f2,Nx*Ny,1)];
else
    if BC_choice==1
    f1 = zeros((Nx+2)*(Ny+1),1);
    for j=0:Ny
       for i=0:Nx+1
           if i==0 || i==Nx+1
               f1((Nx+2)*j+i+1) = k*m_func(Xu(i+1,j+1),Yu(i+1,j+1));
           else
               f1((Nx+2)*j+i+1) = k*dmdx_func(Xu(i+1,j+1),Yu(i+1,j+1));
           end
       end
    end
    elseif BC_choice==2
        f1 = k*reshape(dmdx_func(Xu,Yu),(Nx+1)*(Ny+1),1);  
    end
    f2 = zeros((Nx+1)*(Ny+2),1);
    for j=0:Ny+1
       for i=0:Nx
           if j==0 || j==Ny+1
               f2((Nx+1)*j+i+1) = k*m_func(Xv(i+1,j+1),Yv(i+1,j+1));
           else
               f2((Nx+1)*j+i+1) = k*dmdy_func(Xv(i+1,j+1),Yv(i+1,j+1));
           end
       end
    end
    f = [f1; f2]; 
end

end