function [A] = build_A(Nx,Ny,hx,hy,eta1,eta2,zeta,BC_choice,BC_order)
% author: Mariya Savinov

% Builds the discretization matrix required to solve the force balance
% equation for a compressible Stokesian fluid with 1st and 2nd viscosity
% terms, a contractile element, and drag with an underlying fluid/surface:
%           div(stress) = eta2 * grad(div(u))+eta*L(u) - zeta*u = f
% using finite differences, with a rectangular domain that is either:
%                                   periodic in x and y
%                                   periodic in x, no stress on y boundaries
%                                   no stress at domain boundaries


% Inputs:
%     Nx            : number of x points in the grid
%     Ny            : number of y points in the grid
%     hx            : spatial step size in x
%     hy            : spatial step size in y
%     eta1          : first viscosity coefficient
%     eta2          : second viscosity coefficient
%     zeta          : drag coefficient

% Outputs:
%     A             : matrix A which discretizes eta2 * grad(div(u))+eta*L(u) - zeta*u
    
if BC_choice==0
    %%%%% forming matrix A, which is a discretization of the LHS
    Iy = speye(Ny);  Ix = speye(Nx);  ex = ones(Nx,1); ey = ones(Ny,1);

    % form matrix A1
    T1 = spdiags([ex*(eta1+eta2)/hx^2 (-2*(eta1+eta2)/hx^2 - 2*eta1/hy^2 - zeta)*ex ex*(eta1+eta2)/hx^2],[-1 0 1],Nx,Nx);
    T1(1,end) = (eta1+eta2)/hx^2;               % periodic correction
    T1(end,1) = (eta1+eta2)/hx^2;               % periodic correction

    S = spdiags([ey ey],[-1 1],Ny,Ny);
    A1 = (kron(Iy,T1) + kron(S,Ix*eta1/hy^2));
    A1(end-Nx+1:end,1:Nx) = Ix*eta1/hy^2;       % periodic correction
    A1(1:Nx,end-Nx+1:end) = Ix*eta1/hy^2;       % periodic correction

    % form matrix A2 & A3
    T2 = spdiags([-ex ex],[-1 0],Nx,Nx);
    T2(1,end) = -1;                             % periodic correction

    S = spdiags([-ey ey],[0 1],Ny,Ny);
    A2 = kron(S,T2);
    A2(end-Nx+1:end,1:Nx) =  T2;                % periodic correction

    A2 = eta2/(hx*hy)*A2;
    %%%%%%%%%%%%%%%%%
    T3 = spdiags([-ex ex],[0 1],Nx,Nx);
    T3(end,1) =  1;                             % periodic correction

    S = spdiags([-ey ey],[-1 0],Ny,Ny);
    A3 = kron(S,T3);
    A3(1:Nx,end-Nx+1:end) = -T3;                % periodic correction

    A3 = eta2/(hx*hy)*A3;

    % form matrix A4
    T4 = spdiags([ex*eta1/hx^2 (-2*(eta1+eta2)/hy^2 - 2*eta1/hx^2 - zeta)*ex ex*eta1/hx^2],[-1 0 1],Nx,Nx);
    T4(1,end) = eta1/hx^2;                       % periodic correction
    T4(end,1) = eta1/hx^2;                       % periodic correction

    S = spdiags([ey ey],[-1 1],Ny,Ny);
    A4 = (kron(Iy,T4) + kron(S,Ix*(eta1+eta2)/hy^2));
    A4(end-Nx+1:end,1:Nx) = Ix*(eta1+eta2)/hy^2; % periodic correction
    A4(1:Nx,end-Nx+1:end) = Ix*(eta1+eta2)/hy^2; % periodic correction

elseif BC_choice==1

    %%%% forming matrix A1 %%%%%
    Iy = speye(Ny+1); Ix = speye(Nx+2); ex = ones(Nx+2,1); ey = ones(Ny+1,1);

    % all x-coords, interior y-coords for u
    T1 = spdiags([ex*(eta1+eta2)/hx^2 (-2*(eta1+eta2)/hx^2 - 2*eta1/hy^2 - zeta)*ex ex*(eta1+eta2)/hx^2],[-1 0 1],Nx+2,Nx+2);
    
    % x=a and x=b boundary conditions
    if BC_order==1
        % x=a and x=b boundary conditions
        T1(1,1) = -(eta1+eta2)/hx - hx/(hy^2)*(eta2-eta1);
        T1(1,2) = +(eta1+eta2)/hx;
        T1(end, end-1) = -(eta1+eta2)/hx;
        T1(end,end) = (eta1+eta2)/hx + hx/(hy^2)*(eta2-eta1); 
    else % 2nd order sides
        T1(1,1) = -3*(eta1+eta2)/(2*hx) - hx/(hy^2)*(eta2-eta1);
        T1(1,2) = +2*(eta1+eta2)/hx;
        T1(1,3) = -1/2*(eta1+eta2)/hx;
        T1(end,end-2) = +1/2*(eta1+eta2)/hx;
        T1(end, end-1) = -2*(eta1+eta2)/hx;
        T1(end,end) = +3*(eta1+eta2)/(2*hx) + hx/(hy^2)*(eta2-eta1);
    end
    
    S = spdiags([ey ey],[-1 1],Ny+1,Ny+1);
    I1 = Ix*eta1/hy^2;
    % x=a and x=b boundary conditions
    I1(1,1) = hx/(2*hy^2)*(eta2-eta1);
    I1(end,end) = -hx/(2*hy^2)*(eta2-eta1);
    A1 = (kron(Iy,T1) + kron(S,I1));

    % all x-coords, boundary y-coords for u
    M1 = spdiags([ex*(eta1+eta2)/hx^2 (-2*(eta1+eta2)/hx^2 - eta1/hy^2 - zeta)*ex ex*(eta1+eta2)/hx^2],[-1 0 1],Nx+2,Nx+2);
    
    if BC_order~=2
        M1(1,1) = -(eta1+eta2)/hx;
        M1(1,2) = (eta1+eta2)/hx;
        M1(end,end-1) = -(eta1+eta2)/hx;
        M1(end,end) = (eta1+eta2)/hx;
    else  % corners 2nd order
        M1(1,1) = -3*(eta1+eta2)/(2*hx);
        M1(1,2) = 2*(eta1+eta2)/hx;
        M1(1,3) = -1/2*(eta1+eta2)/hx;
        M1(end,end-2) = 1/2*(eta1+eta2)/hx;
        M1(end,end-1) = -2*(eta1+eta2)/hx;
        M1(end,end) = 3*(eta1+eta2)/(2*hx);
    end
    
    M2 = I1;
    M2(1,1) = 0;
    M2(end,end) = 0;

    A1(1:Nx+2,1:Nx+2) = M1;
    A1(1:Nx+2,Nx+3:2*Nx+4) = M2;
    A1(end-Nx-1:end,end-Nx-1:end) = M1;
    A1(end-Nx-1:end,end-2*Nx-3:end-Nx-2) = M2;


    %%%% forming matrix A2 %%%%%
    % all x-coords, interior y-coords for u 
    T2 = spdiags([-ex ex],[-1 0],Nx+2,Nx+1);
    T2 = (eta2)/hx/hy*T2;
    % x=a and x=b boundary conditions
    T2(1,1) = (eta2-eta1)/hy;
    T2(end,end) = (eta2-eta1)/hy;

    S = spdiags([-ey ey],[0 1],Ny+1,Ny+2);
    A2 = kron(S,T2);

    % all x-coords, boundary y-coords for u
    M3 = spdiags([-ex +ex],[-1 0],Nx+2,Nx+1);
    M4 = (eta2)/(hx*hy)*M3;
    M3 = (eta1-eta2)/(hx*hy)*M3;

    M3(1,1) = -3/(2*hy)*(eta2-eta1);
    M3(1,2) = (eta2-eta1)/(2*hy);
    M3(end,end) = -3/(2*hy)*(eta2-eta1);
    M3(end,end-1) = (eta2-eta1)/(2*hy);
    M4(1,1) = +3/(2*hy)*(eta2-eta1);
    M4(1,2) = -(eta2-eta1)/(2*hy);
    M4(end,end) = +3/(2*hy)*(eta2-eta1);
    M4(end,end-1) = -(eta2-eta1)/(2*hy);

    A2(1:Nx+2,1:Nx+1) = M3;
    A2(1:Nx+2,Nx+2:2*Nx+2) = M4;
    A2(end-Nx-1:end,end-Nx:end) = -M3;
    A2(end-Nx-1:end,end-2*Nx-1:end-Nx-1) = -M4;

    %%%% forming matrix A4 %%%%%
    Iy = speye(Ny+2); Ix = speye(Nx+1); ex = ones(Nx+1,1); ey = ones(Ny+2,1);

    T4 = spdiags([ex*(eta1/hx^2) (-2*(eta1+eta2)/(hy^2) - 2*eta1/(hx^2) - zeta)*ex ex*(eta1/hx^2)],[-1 0 1],Nx+1,Nx+1);
    % adjustments for ghost points near x=a and x=b
    T4(1,1) = -2*(eta1+eta2)/(hy^2)-eta1/hx^2-zeta;
    T4(end,end) = -2*(eta1+eta2)/(hy^2)-eta1/hx^2-zeta;


    S = spdiags([ey ey],[-1 1],Ny+2,Ny+2);
    I4 = Ix*(eta1+eta2)/hy^2; 

    A4 = (kron(Iy,T4) + kron(S,I4));

    % y=c and y=d boundary conditions
    
    if BC_order==1
        M5 = spdiags([ex*(hy/(2*hx^2)*(eta2-eta1)) ex*(-(eta2+eta1)/hy-(eta2-eta1)*hy/(hx^2)) ex*(hy/(2*hx^2)*(eta2-eta1))], [-1 0 1], Nx+1,Nx+1);
        M5(1,1) = -(eta2+eta1)/hy;
        M5(1,2) = 0;
        M5(end,end) = -(eta2+eta1)/hy;
        M5(end,end-1) = 0;

        M6 = Ix*(eta2+eta1)/hy;

        A4(1:Nx+1,1:Nx+1) = M5;
        A4(1:Nx+1,Nx+2:2*Nx+2) = M6;
        A4(end-Nx:end,end-Nx:end) = -M5;
        A4(end-Nx:end,end-2*Nx-1:end-Nx-1) = -M6;
    elseif BC_order==21
        M5 = spdiags([ex*(hy/(2*hx^2)*(eta2-eta1)) ex*(-3*(eta2+eta1)/(2*hy)-(eta2-eta1)*hy/(hx^2)) ex*(hy/(2*hx^2)*(eta2-eta1))], [-1 0 1], Nx+1,Nx+1);
        M5(1,1) = -(eta2+eta1)/hy;
        M5(1,2) = 0;
        M5(end,end) = -(eta2+eta1)/hy;
        M5(end,end-1) = 0;

        M6 = Ix*2*(eta2+eta1)/hy;
        M7 = Ix*(-1/2)*(eta2+eta1)/hy;
        M6(1,1) = (eta2+eta1)/hy;
        M6(end,end) = (eta2+eta1)/hy;
        M7(1,1)=0;
        M7(end,end)=0;
        
        A4(1:Nx+1,1:Nx+1) = M5;
        A4(1:Nx+1,Nx+2:2*Nx+2) = M6;
        A4(1:Nx+1,2*Nx+3:3*Nx+3) = M7;
        A4(end-Nx:end,end-Nx:end) = -M5;
        A4(end-Nx:end,end-2*Nx-1:end-Nx-1) = -M6;
        A4(end-Nx:end,end-3*Nx-2:end-2*Nx-2) = -M7;
    else % all 2nd order
        M5 = spdiags([ex*(hy/(2*hx^2)*(eta2-eta1)) ex*(-3*(eta2+eta1)/(2*hy)-(eta2-eta1)*hy/(hx^2)) ex*(hy/(2*hx^2)*(eta2-eta1))], [-1 0 1], Nx+1,Nx+1);
        M5(1,1) = -3*(eta2+eta1)/(2*hy);
        M5(1,2) = 0;
        M5(end,end) = -3*(eta2+eta1)/(2*hy);
        M5(end,end-1) = 0;

        M6 = Ix*2*(eta2+eta1)/hy;
        M7 = Ix*(-1/2)*(eta2+eta1)/hy;

        A4(1:Nx+1,1:Nx+1) = M5;
        A4(1:Nx+1,Nx+2:2*Nx+2) = M6;
        A4(1:Nx+1,2*Nx+3:3*Nx+3) = M7;
        A4(end-Nx:end,end-Nx:end) = -M5;
        A4(end-Nx:end,end-2*Nx-1:end-Nx-1) = -M6;
        A4(end-Nx:end,end-3*Nx-2:end-2*Nx-2) = -M7;
    end
    
    %%%%%%% forming matrix A3

    % all x-coords, interior y-coords for u 
    T3 = spdiags([-ex ex],[0 1],Nx+1,Nx+2);
    T3 = (eta2)/hx/hy*T3;
    % x=a and x=b ghost points  (checked)
    T3(1,1) = +(eta1-eta2)/(hx*hy);
    T3(1,2) = + eta2/(hx*hy);    
    T3(end,end) = - (eta1-eta2)/(hx*hy); 
    T3(end,end-1) = - eta2/(hx*hy);     

    S = spdiags([-ey ey],[-1 0],Ny+2,Ny+1);
    A3 = kron(S,T3);

    % y=c and y=d boundary conditions
    M9 = spdiags([-ex ex],[0 1],Nx+1,Nx+2);
    M9 = (eta2-eta1)/(hx)*M9;
    M9(1,1) = -3/(2*hx)*(eta2-eta1);
    M9(1,2) = +3/(2*hx)*(eta2-eta1);
    M9(end,end) = +3/(2*hx)*(eta2-eta1);
    M9(end,end-1) = -3/(2*hx)*(eta2-eta1);

    M10 = sparse(Nx+1,Nx+2);
    M10(1,1) = 1/(2*hx)*(eta2-eta1);
    M10(1,2) = -1/(2*hx)*(eta2-eta1);
    M10(end,end) = -1/(2*hx)*(eta2-eta1);
    M10(end,end-1) = 1/(2*hx)*(eta2-eta1);


    A3(1:Nx+1,1:Nx+2) = M9;
    A3(1:Nx+1,Nx+3:2*Nx+4) = M10;
    A3(end-Nx:end,end-Nx-1:end) = M9;
    A3(end-Nx:end,end-2*Nx-3:end-Nx-2) = M10;

elseif BC_choice==2
    %%%% forming matrix A1 %%%%%
    Iy = speye(Ny+1); Ix = speye(Nx+1); ex = ones(Nx+1,1); ey = ones(Ny+1,1);

    % all x-coords, interior y-coords for u
    T1 = spdiags([ex*(eta1+eta2)/hx^2 (-2*(eta1+eta2)/hx^2 - 2*eta1/hy^2 - zeta)*ex ex*(eta1+eta2)/hx^2],[-1 0 1],Nx+1,Nx+1);
    % periodic in x conditions
    T1(1,end) = (eta1+eta2)/hx^2;
    T1(end,1) = (eta1+eta2)/hx^2;

    S = spdiags([ey ey],[-1 1],Ny+1,Ny+1);
    I1 = Ix*eta1/hy^2;

    A1 = (kron(Iy,T1) + kron(S,I1));

    % all x-coords, boundary y-coords for u
    M1 = spdiags([ex*(eta1+eta2)/hx^2 (-2*(eta1+eta2)/hx^2 - eta1/hy^2 - zeta)*ex ex*(eta1+eta2)/hx^2],[-1 0 1],Nx+1,Nx+1);
    M1(1,end) = (eta1+eta2)/hx^2;
    M1(end,1) = (eta1+eta2)/hx^2;

    M2 = I1;

    A1(1:Nx+1,1:Nx+1) = M1;
    A1(1:Nx+1,Nx+2:2*Nx+2) = M2;
    A1(end-Nx:end,end-Nx:end) = M1;
    A1(end-Nx:end,end-2*Nx-1:end-Nx-1) = M2;


    %%%% forming matrix A2 %%%%%
    % all x-coords, interior y-coords for u 
    T2 = spdiags([-ex ex],[-1 0],Nx+1,Nx+1);
    % periodic in x
    T2(1,end) = -1;

    T2 = (eta2)/hx/hy*T2;

    S = spdiags([-ey ey],[0 1],Ny+1,Ny+2);
    A2 = kron(S,T2);

    % all x-coords, boundary y-coords for u
    M3 = T2*(eta1-eta2)/(eta2);

    A2(1:Nx+1,1:Nx+1) = M3;
    A2(end-Nx:end,end-Nx:end) = -M3;

    %%%% forming matrix A4 %%%%%
    Iy = speye(Ny+2); Ix = speye(Nx+1); ex = ones(Nx+1,1); ey = ones(Ny+2,1);

    T4 = spdiags([ex*(eta1/hx^2) (-2*(eta1+eta2)/(hy^2) - 2*eta1/(hx^2) - zeta)*ex ex*(eta1/hx^2)],[-1 0 1],Nx+1,Nx+1);
    % periodic in x
    T4(1,end) = (eta1/hx^2);
    T4(end,1) = (eta1/hx^2);

    S = spdiags([ey ey],[-1 1],Ny+2,Ny+2);
    I4 = Ix*(eta1+eta2)/hy^2; 

    A4 = (kron(Iy,T4) + kron(S,I4));


    if BC_order==1 % first order
        % y=c and y=d boundary conditions
        M5 = spdiags([ex*(hy/(2*hx^2)*(eta2-eta1)) ex*(-(eta2+eta1)/hy-(eta2-eta1)*hy/(hx^2)) ex*(hy/(2*hx^2)*(eta2-eta1))], [-1 0 1], Nx+1,Nx+1);
        M5(1,end) = (hy/(2*hx^2)*(eta2-eta1));
        M5(end,1) = (hy/(2*hx^2)*(eta2-eta1));

        M6 = Ix*(eta2+eta1)/hy;

        A4(1:Nx+1,1:Nx+1) = M5;
        A4(1:Nx+1,Nx+2:2*Nx+2) = M6;
        A4(end-Nx:end,end-Nx:end) = -M5;
        A4(end-Nx:end,end-2*Nx-1:end-Nx-1) = -M6;
    elseif BC_order==2
        % y=c and y=d boundary conditions
        M5 = spdiags([ex*(hy/(2*hx^2)*(eta2-eta1)) ex*(-3*(eta2+eta1)/(2*hy)-(eta2-eta1)*hy/(hx^2)) ex*(hy/(2*hx^2)*(eta2-eta1))], [-1 0 1], Nx+1,Nx+1);
        M5(1,end) = (hy/(2*hx^2)*(eta2-eta1));
        M5(end,1) = (hy/(2*hx^2)*(eta2-eta1));

        M6 = Ix*2*(eta2+eta1)/hy;

        M7 = Ix*(-1/2)*(eta2+eta1)/hy;

        A4(1:Nx+1,1:Nx+1) = M5;
        A4(1:Nx+1,Nx+2:2*Nx+2) = M6;
        A4(1:Nx+1,2*Nx+3:3*Nx+3) = M7;
        A4(end-Nx:end,end-Nx:end) = -M5;
        A4(end-Nx:end,end-2*Nx-1:end-Nx-1) = -M6;
        A4(end-Nx:end,end-3*Nx-2:end-2*Nx-2) = -M7;
    end

    %%%%%%% forming matrix A3

    % all x-coords, interior y-coords for u 
    T3 = spdiags([-ex ex],[0 1],Nx+1,Nx+1);
    % periodic in x
    T3(end,1) = 1;
    T3 = (eta2)/hx/hy*T3;

    S = spdiags([-ey ey],[-1 0],Ny+2,Ny+1);
    A3 = kron(S,T3);

    % y=c and y=d boundary conditions
    M9 = T3*hy*(eta2-eta1)/eta2;

    M10 = sparse(Nx+1,Nx+1);

    A3(1:Nx+1,1:Nx+1) = M9;
    A3(1:Nx+1,Nx+2:2*Nx+2) = M10;
    A3(end-Nx:end,end-Nx:end) = M9;
    A3(end-Nx:end,end-2*Nx-1:end-Nx-1) = M10;
end


% form matrix A from blocks A1,A23,A4
A = [A1 A2;...
     A3 A4];

