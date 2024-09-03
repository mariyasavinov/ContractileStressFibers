function [Fibers_array,u_array,v_array,m_array,b_array,t_array,iterK0] = fiber_2D_solver(XYuvm,hx,hy,Fibers,nsteps,dt,eta1,eta2,zeta,d0,alpha,beta,gamma,f0,m0,b0,global_vars,Fiber_vars)
% author: Mariya Savinov

% Solves the force balance equation for a compressible Stokesian fluid
% with 1st and 2nd viscosity terms, a contractile element, and drag with
% underlying fluid/surface:
%           div(stress) = eta2 * grad(div(u))+eta*L(u) - zeta*u = f
% on 2D periodic rectangular domain, using finite differences


% Inputs:
%     XYuv          : x,y values from Eulerian meshgrid, such that X(i,j),Y(i,j) gives the x_i, y_j point
%     hx            : spatial step size in x
%     hy            : spatial step size in y
%     eta1          : first viscosity coefficient
%     eta2          : second viscosity coefficient
%     zeta          : drag coefficient
%     f             : source/forcing on RHS (array of all elements)

% Outputs:
%     u      :  solution on [X,Y] for x-direction of velocity field
%     v      :  solution on [X,Y] for y-direction of velocity field


ActinFluid = global_vars(1); BC_choice = global_vars(2); BC_order = global_vars(3); PlotWhileRun = global_vars(4);
FixedPointIteration = global_vars(5); FPItol = global_vars(6); AdvFlux_choice = global_vars(7); FreeMyosin = global_vars(8);

FiberLoop = Fiber_vars(1); FiberContractile = Fiber_vars(2); FiberAdhesion = Fiber_vars(3);
FiberViscous = Fiber_vars(4); FiberSlip = Fiber_vars(5); FiberBending = Fiber_vars(6);
LaserAblation = Fiber_vars(7); Indx_Laser = Fiber_vars(8); T_Laser = Fiber_vars(9);

if ActinFluid==0
    gamma = gamma*Fibers{1}{2}(1);
end
 

X = XYuvm{1}; Y = XYuvm{2};
Xu = XYuvm{3}; Yu = XYuvm{4};
Xv = XYuvm{5}; Yv = XYuvm{6};
Xm = XYuvm{7}; Ym = XYuvm{8};

dtsave = 0.01;
dtsave_indx = round(dtsave/dt);

% Get bounds of domain
La = min(X,[],'all'); Lb = max(X,[],'all')+hx; Lc = min(Y,[],'all'); Ld = max(Y,[],'all')+hy;

% Get number of X pts and number of Y pts
if BC_choice==0
    Nx = size(X,1);
    Ny = size(X,2);
    mNx = Nx; mNy = Ny;
else
   Nx = size(X,1)-2;
   Ny = size(X,2)-2;
   mNx = Nx+1; mNy = Ny+1;
end
% 
% Nb = size(XYfiber,1);

% build and save A for computations involving the actin bulk fluid
[A] = build_A(Nx,Ny,hx,hy,eta1,eta2,zeta,BC_choice,BC_order);

% build and save L for computations if you have diffusing free myosin
if FreeMyosin==1
Ix = speye(mNx); Iy = speye(mNy); ex = ones(mNx,1); ey = ones(mNy,1);
T = spdiags([ex/hx^2 -2*(1/hx^2+1/hy^2)*ex ex/hx^2],[-1 0 1],mNx,mNx);
if BC_choice==0     % periodic
    T(1,end) = 1/hx^2;
    T(end,1) = 1/hx^2;
elseif BC_choice==1 % no flux
    T(1,1) = T(1,1) + 1/hx^2;
    T(end,end) = T(end,end) + 1/hx^2;
end
S = spdiags([ey ey],[-1 1],mNy,mNy);
L = (kron(Iy,T) + kron(S,Ix/hy^2));
if BC_choice==0     % periodic
    L(end-mNx+1:end,1:mNx) = Ix/hy^2;
    L(1:mNx,end-mNx+1:end) = Ix/hy^2;
elseif BC_choice==1 % no flux
    T = T + diag((1/hy^2)*ex);
    L(1:mNx,1:mNx) = T;
    L(end-mNx+1:end,end-mNx+1:end) = T;
end
L = d0*L;
I = speye(mNx*mNy);
end

% tracking arrays
m_array = {}; b_array = {};
% XYfiber_array = {};
% Ufiber_array = {}; Vfiber_array = {};
Fibers_array = {};
u_array = {}; v_array = {};

t_array = [];


% initial condition
t = 0; m = m0; b = b0;

if PlotWhileRun==1
if FreeMyosin~=1
figure('position', [100, 100, 600, 600]);
else
figure('position', [100, 100, 1600, 400]);
subaxis(1,2,1)
end
hold on
hfig = pcolor(Xm,Ym,m);
set(hfig,'EdgeColor','none')
colorbar
for k = 1:length(Fibers)
    XYfiber = Fibers{k}{1};
    Zfiber = Fibers{k}{4};
    if Fibers{k}{8} == 1
    plot([XYfiber(:,1);XYfiber(1,1)],[XYfiber(:,2);XYfiber(1,2)],'w.-','LineWidth',3)
    else
    plot(XYfiber(:,1),XYfiber(:,2),'w.-','LineWidth',3)
    end
    if isempty(Zfiber) ~= 1
        plot(Zfiber(:,1),Zfiber(:,2),'r*','MarkerSize',3)
    end
end
axis equal
if FreeMyosin==1
subaxis(1,2,2)
hfig = pcolor(Xm,Ym,b);
set(hfig,'EdgeColor','none')
colorbar
ylabel('y','Fontsize',14)
xlabel('x','Fontsize',14)
clim([0,1])
axis([La Lb Lc Ld])
axis equal


end
arrow_num = 22;
x_spaced = La:(Lb-La)/arrow_num:Lb;
y_spaced = Lc:(Ld-Lc)/arrow_num:Ld;
[Xq,Yq] = meshgrid(x_spaced,y_spaced);
end

wbar = waitbar(0,'Please wait...');
pause(.05)

% initial guess for fiber velocity, to be used for forces that require it
for k=1:length(Fibers)
    XYfiber = Fibers{k}{1};
    Ufiber = zeros(size(XYfiber,1),1);
    Vfiber = zeros(size(XYfiber,1),1);
    UVfiber = [Ufiber, Vfiber];
    Fibers{k}{3} = UVfiber;
end

iterK0 = 0; % how many iterations first step took
stopSIM = 0;

for n=1:nsteps+1   % +1 to make sure you get configuration AT final time also
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laser Ablation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if LaserAblation==1 && (abs(t-T_Laser)<dt/10 || (t<T_Laser && t+dt>T_Laser))

    % cut first fiber
    XYfiber = Fibers{1}{1}; Nb = size(XYfiber,1);
    hq_array = Fibers{1}{2};
    UVfiber = Fibers{1}{3}; Ufiber = UVfiber(:,1); Vfiber = UVfiber(:,2);
    Zfiber = Fibers{1}{4};
    fiberparams = Fibers{1}{7};
    kBfiber = fiberparams{2};

    % Second Half
    XYfiber = XYfiber(Indx_Laser+1:end,:);
    Ufiber = Ufiber(Indx_Laser+1:end); Vfiber = Vfiber(Indx_Laser+1:end);
    UVfiber = [Ufiber, Vfiber];
    hq_array = hq_array(Indx_Laser+1:end);
    Nb = Nb-Indx_Laser;

    % rebuild and save B for bending forces
    if FiberBending==1
    [B] = build_B(Nb,hq_array(1),kBfiber,FiberLoop);
    if FiberViscous==1 || FiberContractile==3 %%% if you are doing FPI
        % reorder so it is X1;Y1;X2;Y2;...
        B_reorder = sparse(2*Nb,2*Nb);
        B_reorder(1:2:end,1:2:end) = B(1:Nb,1:Nb);
        B_reorder(2:2:end,2:2:end) = B(Nb+1:end,Nb+1:end);
        B = B_reorder;
    end
    else
    B = sparse(2*Nb,2*Nb);
    end
    % C and D matrices for Slip with FPI
    C = sparse(2*Nb,2*Nb); 
    D = sparse(2*Nb,2*Nb);

    % identity matrix
    I2Nb = speye(2*Nb);

    NewIndx = length(Fibers);
    Fibers{NewIndx+1}{1} = XYfiber;
    Fibers{NewIndx+1}{2} = hq_array;
    Fibers{NewIndx+1}{3} = UVfiber;
    if isempty(Zfiber)~=1
        ZfiberNew = Zfiber(find(Zfiber(:,3)>Indx_Laser),:);
        ZfiberNew(:,3) = ZfiberNew(:,3) - Indx_Laser;
    Fibers{NewIndx+1}{4} = ZfiberNew;
    % Fibers{NewIndx+1}{4} = [Zfiber(2,1),Zfiber(2,2),Zfiber(2,3)-Indx_Laser];
    else; Fibers{NewIndx+1}{4} = []; end
    Fibers{NewIndx+1}{5} = {B,C,D};
    Fibers{NewIndx+1}{6} = I2Nb;
    Fibers{NewIndx+1}{7} = fiberparams;
    Fibers{NewIndx+1}{8} = 0; % not a loop
    Fibers{NewIndx+1}{9} = NaN;

%%%%%%%%%%%%%%%%%%%%%
    XYfiber = Fibers{1}{1}; hq_array = Fibers{1}{2};
    UVfiber = Fibers{1}{3}; Ufiber = UVfiber(:,1); Vfiber = UVfiber(:,2);

    % First Half
    XYfiber = XYfiber(1:Indx_Laser,:);
    Ufiber = Ufiber(1:Indx_Laser); Vfiber = Vfiber(1:Indx_Laser);
    UVfiber = [Ufiber,Vfiber];
    hq_array = hq_array(1:Indx_Laser-1);
    Nb = Indx_Laser;

    % rebuild and save B for bending forces
    if FiberBending==1
    [B] = build_B(Nb,hq_array(1),kBfiber,FiberLoop);
    if FiberViscous==1  %%% if you are doing FPI
        % reorder so it is X1;Y1;X2;Y2;...
        B_reorder = sparse(2*Nb,2*Nb);
        B_reorder(1:2:end,1:2:end) = B(1:Nb,1:Nb);
        B_reorder(2:2:end,2:2:end) = B(Nb+1:end,Nb+1:end);
        B = B_reorder;
    end
    else
    B = sparse(2*Nb,2*Nb);
    end
    % C and D matrices for Slip with FPI
    C = sparse(2*Nb,2*Nb); 
    D = sparse(2*Nb,2*Nb);

    % identity matrix
    I2Nb = speye(2*Nb);

    NewIndx = length(Fibers);
    Fibers{1}{1} = XYfiber;
    Fibers{1}{2} = hq_array;
    Fibers{1}{3} = UVfiber;
    if isempty(Zfiber)~=1
        ZfiberNew = Zfiber(find(Zfiber(:,3)<=Indx_Laser),:);
    Fibers{1}{4} = ZfiberNew;
    % Fibers{1}{4} = [Zfiber(1,:)];
    end
    Fibers{1}{5} = {B,C,D};
    Fibers{1}{6} = I2Nb;
    % Fibers{NewIndx+1}{7} = fiberparams;
    % Fibers{NewIndx+1}{8} = 0; % not a loop
    % Fibers{NewIndx+1}{9} = NaN;
    LaserAblation = 0;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% taking 1 timestep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = t+dt;

if FixedPointIteration==1 % if FiberViscous==1 or FiberContractile==3

% build matrices C and D for FPI
for k=1:length(Fibers)
    XYfiber = Fibers{k}{1}; Nb = size(XYfiber,1);
    hq_array = Fibers{k}{2};
    xi = Fibers{k}{7}{4};
    contrConst = Fibers{k}{7}{5};
    BCD = Fibers{k}{5};
    if FiberSlip==1
        if FiberViscous==1
            [Av] = build_AFv(XYfiber,hq_array);
            C = xi*Av;
            BCD{2} = C;
        end
        if FiberContractile==3
            % Cm = build_AFv(XYfiber,ones(size(hq_array)));
            Cm = build_AFv(XYfiber,hq_array);
            Fs = contrConst(3); v0 = contrConst(4);
            D = Fs/v0*Cm;
            BCD{3} = D;
        end
    end
    Fibers{k}{5} = BCD;
end

maxIter = 150;

% save old XYfiber before updating
for k=1:length(Fibers)
    XYfiberOld = Fibers{k}{1};
    Fibers{k}{10}=XYfiberOld;
end

iterBreak = 0;
for iterK = 1:maxIter

    Ffiber_spread_All = [];

    for k=1:length(Fibers)
        XYfiber = Fibers{k}{1};
        XYfiberOld = Fibers{k}{10};
        hq_array = Fibers{k}{2};
        UVfiber = Fibers{k}{3}; Ufiber = UVfiber(:,1); Vfiber = UVfiber(:,2);
        Zfiber = Fibers{k}{4};
        fiberparams = Fibers{k}{7};
        kfiber = fiberparams{1}; kBfiber = fiberparams{2}; kadh = fiberparams{3};
        xi = fiberparams{4}; contrConst = fiberparams{5};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % given configuration, solve for fluid velocity
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Ffiber_spread, felastic, fbending, fviscous, fadhesion, fcontractile]= get_Ffiber_spread(XYuvm,XYfiberOld,Ufiber,Vfiber,Zfiber,hq_array,kfiber,kBfiber,xi,contrConst,kadh,BC_choice,Fiber_vars);
        
        if isempty(Ffiber_spread_All); Ffiber_spread_All = Ffiber_spread;
        else; Ffiber_spread_All = Ffiber_spread_All + Ffiber_spread; end

        % temporarily save felastic, fadhesion, and fcontractile for timestepping
        Fibers{k}{9} = {felastic, fadhesion, fcontractile, fbending, fviscous};
    
        UVfiberOld = UVfiber; Fibers{k}{11} = UVfiberOld;
        XYfiber_iterOld = XYfiber; Fibers{k}{12} = XYfiber_iterOld;
    
    end


    % if immersed in fluid
    if ActinFluid==1
    % get force due to myosin
    fmyosin = myosinForce_discrete(XYuvm,1,m,hx,hy,BC_choice);
    % solve fluid equation
    [u,v] = fluid_2D_solver(XYuvm,A,f0,fmyosin,Ffiber_spread_All,BC_choice);
    end

    for k = 1:length(Fibers)
        XYfiber = Fibers{k}{1}; Nb = size(XYfiber,1);
        XYfiberOld = Fibers{k}{10};
        BCD = Fibers{k}{5}; I2Nb = Fibers{k}{6};
        B = BCD{1}; C = BCD{2}; D = BCD{3};
        felastic = Fibers{k}{9}{1};
        fadhesion = Fibers{k}{9}{2};
        fcontractile = Fibers{k}{9}{3};
        UVfiberOld = Fibers{k}{11};
        XYfiber_iterOld = Fibers{k}{12};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get component from fluid, if applicable
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ActinFluid==1
        % interpolate velocity to the fiber
        Ufiber = u_interp1D(Xu,Yu,XYfiberOld,u); % changed XYfiber --> XYfiberOld
        Vfiber = u_interp1D(Xv,Yv,XYfiberOld,v);
        else
            Ufiber = 0;
            Vfiber = 0;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update position of fiber
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if FiberSlip==1

            RHS = zeros(2*Nb,1);
            % RHS(1:2:end) = Ufiber + 1/gamma*felastic(:,1) + 1/gamma*fadhesion(:,1) + 1/gamma*fbending(:,1);
            % RHS(2:2:end) = Vfiber + 1/gamma*felastic(:,2) + 1/gamma*fadhesion(:,2) + 1/gamma*fbending(:,2);
            XYfiber_vec = zeros(2*Nb,1); XYfiber_vec(1:2:end) = XYfiberOld(:,1); XYfiber_vec(2:2:end) = XYfiberOld(:,2);
            RHS(1:2:end) = XYfiberOld(:,1)+dt*(Ufiber + 1/gamma*felastic(:,1) + 1/gamma*fadhesion(:,1));
            RHS(2:2:end) = XYfiberOld(:,2)+dt*(Vfiber + 1/gamma*felastic(:,2) + 1/gamma*fadhesion(:,2));
            % RHS = RHS - 1/gamma*C*XYfiber_vec;
            RHS = RHS - 1/gamma*C*XYfiber_vec - 1/gamma*D*XYfiber_vec;
            if FiberContractile~=3
                RHS(1:2:end) = RHS(1:2:end) + dt/gamma*fcontractile(:,1);
                RHS(2:2:end) = RHS(2:2:end) + dt/gamma*fcontractile(:,2);
            else
                tau = (XYfiberOld(2:end,:)-XYfiberOld(1:end-1,:))./vecnorm(XYfiberOld(2:end,:)-XYfiberOld(1:end-1,:),2,2);
                Z = Fs.*tau; Z = [Z;[0,0]];
                if FiberLoop == 1
                tau_conn = (XYfiberOld(1,:)-XYfiberOld(end,:))./vecnorm(XYfiberOld(1,:)-XYfiberOld(end,:),2,2);
                Z(end,:) = Fs.*tau_conn;
                end
                FcontrEX = (Z - [Z(end,:); Z(1:end-1,:)]); % Fs(tau(k+1/2) - tau(k-1/2)
                RHS(1:2:end) = RHS(1:2:end) + dt/gamma*FcontrEX(:,1);
                RHS(2:2:end) = RHS(2:2:end) + dt/gamma*FcontrEX(:,2);
            end

            LHS = speye(2*Nb) - 1/gamma*C - dt/gamma*B - 1/gamma*D;

            if n==1 && iterK==1
                fprintf('h=%.2e, condLHS = %.2e\n',hq_array(1),cond(LHS))
            end
            xLHS = LHS\RHS;

            XYfiber(:,1) = xLHS(1:2:end); XYfiber(:,2) = xLHS(2:2:end);

            % approximate fiber velocity to track
            Ufiber = (XYfiber(:,1)-XYfiberOld(:,1))/dt; Vfiber=(XYfiber(:,2)-XYfiberOld(:,2))/dt;
            UVfiber = [Ufiber, Vfiber];
        else 
            UVfiber = [Ufiber, Vfiber];
            % Forward Euler step
            XYfiber = XYfiberOld + dt*[Ufiber,Vfiber];
        end

        Fibers{k}{1} = XYfiber;
        Fibers{k}{3} = UVfiber;
    
        if norm(vecnorm(UVfiberOld-UVfiber,2,2))/norm(vecnorm(UVfiber,2,2))<FPItol
            iterBreak = iterBreak+1;
        end
    end

    if iterBreak == length(Fibers)
        waitbar(n/(nsteps),wbar,sprintf('Simulation running, FPI converged at iteration %d\n',n,iterK));
        if n==1 % save first iteration count
            iterK0 = iterK;
        end
        break
    else
        iterBreak = 0;
    end

    if iterK == maxIter
        disp('Exceeded max iterations for fixed point iteration')
        iterK0 = NaN;
        stopSIM = 1;
        break
    end
end

else
waitbar(n/(nsteps),wbar,'Simulation running');

Ffiber_spread_All = [];

for k=1:length(Fibers)

    XYfiber = Fibers{k}{1};
    XYfiberOld = XYfiber; Fibers{k}{10}=XYfiberOld;
    hq_array = Fibers{k}{2};
    UVfiber = Fibers{k}{3}; Ufiber = UVfiber(:,1); Vfiber = UVfiber(:,2);
    Zfiber = Fibers{k}{4};
    fiberparams = Fibers{k}{7};
    kfiber = fiberparams{1}; kBfiber = fiberparams{2}; kadh = fiberparams{3};
    xi = fiberparams{4}; contrConst = fiberparams{5};
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % given configuration, get forces then solve for fluid velocity and update fiber position 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Ffiber_spread, felastic, ~, ~, fadhesion, fcontractile]= get_Ffiber_spread(XYuvm,XYfiber,Ufiber,Vfiber,Zfiber,hq_array,kfiber,kBfiber,xi,contrConst,kadh,BC_choice,Fiber_vars);

    if isempty(Ffiber_spread_All); Ffiber_spread_All = Ffiber_spread;
    else; Ffiber_spread_All = Ffiber_spread_All + Ffiber_spread; end

    % temporarily save felastic, fadhesion, and fcontractile for timestepping
    Fibers{k}{9} = {felastic, fadhesion, fcontractile};
end

% if immersed in fluid
if ActinFluid==1
% get force due to myosin
fmyosin = myosinForce_discrete(XYuvm,1,m,hx,hy,BC_choice);
% solve fluid equation
[u,v] = fluid_2D_solver(XYuvm,A,f0,fmyosin,Ffiber_spread_All,BC_choice);
end

for k = 1:length(Fibers)
    XYfiber = Fibers{k}{1}; Nb = size(XYfiber,1);
    XYfiberOld = Fibers{k}{10};
    B = Fibers{k}{5}{1}; I2Nb = Fibers{k}{6};
    felastic = Fibers{k}{9}{1};
    fadhesion = Fibers{k}{9}{2};
    fcontractile = Fibers{k}{9}{3};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get component from fluid, if applicable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ActinFluid==1
    % interpolate velocity to the fiber
    Ufiber = u_interp1D(Xu,Yu,XYfiber,u);
    Vfiber = u_interp1D(Xv,Yv,XYfiber,v);
    else
    Ufiber = zeros(size(XYfiber,1),1);
    Vfiber = zeros(size(XYfiber,1),1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update fiber position 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if FiberSlip==1 
        % add correction/slip term
        RHS = zeros(2*Nb,1);
        RHS(1:Nb) = XYfiber(:,1)+dt*(Ufiber + 1/gamma*felastic(:,1) + 1/gamma*fcontractile(:,1) + 1/gamma*fadhesion(:,1));
        RHS(Nb+1:end) = XYfiber(:,2) + dt*(Vfiber + 1/gamma*felastic(:,2) + 1/gamma*fcontractile(:,2) + 1/gamma*fadhesion(:,2));
        XYfiber_vec = (I2Nb-dt/gamma*B)\RHS;
        XYfiber(:,1) = XYfiber_vec(1:Nb); XYfiber(:,2) = XYfiber_vec(Nb+1:end);
        % approximate fiber velocity to track
        Ufiber = (XYfiber(:,1)-XYfiberOld(:,1))/dt; Vfiber=(XYfiber(:,2)-XYfiberOld(:,2))/dt;
    else
        % Forward Euler step
        XYfiber = XYfiber + dt*[Ufiber,Vfiber];
        % note: if FiberSlip=0 and ActinFluid=0, fiber will not move
    end
    
    Fibers{k}{1} = XYfiber;
    Fibers{k}{3} = [Ufiber,Vfiber];

end
end

% save fluid velocity, myosin, current fiber configuration and velocity
if mod(n-1,dtsave_indx)==0 || n==nsteps+1
    m_array{end+1} = m; b_array{end+1} = b;
    if ActinFluid==1;   u_array{end+1} = u; v_array{end+1} = v; end
    
    % Ufiber_array{end+1} = Ufiber; Vfiber_array{end+1} = Vfiber;
    % XYfiber_array{end+1} = XYfiberOld;
    
    Tindx = length(Fibers_array)+1;

    for k = 1:length(Fibers)
        % Fibers_XYfiber_array{Tindx}{k} = Fibers{k}{10}; %XYfiberOld
        % Fibers_UVfiber_array{Tindx}{k} = Fibers{k}{3};
        % Fibers_hq_array{Tindx}{k} = Fibers{k}{2};
        Fibers_array{Tindx}{k}{1} = Fibers{k}{10}; % XYfiberOld;
        Fibers_array{Tindx}{k}{2} = Fibers{k}{2}; % hq_array
        Fibers_array{Tindx}{k}{3} = Fibers{k}{3}; % UVfiber
        Fibers_array{Tindx}{k}{4} = Fibers{k}{4}; % Zfiber
        Fibers_array{Tindx}{k}{5} = Fibers{k}{5}; % fiberparams
    end

    t_array(end+1) = t-dt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update myosin density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ActinFluid==1
if FreeMyosin~=1
m = myosin_2D_TimeStep(u,v,m,hx,hy,dt,zeros(size(m)),AdvFlux_choice,BC_choice);
else
m = myosin_2D_TimeStep(u,v,m,hx,hy,dt,-alpha*m+beta*b,AdvFlux_choice,BC_choice);

RHS = b + dt*alpha*m - dt*beta*b;
RHS_vec = reshape(RHS,mNx*mNy,1) + dt/2*L*reshape(b,mNx*mNy,1);
b_vec = (I-dt/2*L)\RHS_vec;
b = reshape(b_vec,mNx,mNy);
end
end

if PlotWhileRun==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (mod(t,0.001)>=0 && mod(t,0.001)<=dt && t>0.001+2*dt) 

if ActinFluid==1
u_spaced = interp2(Yu,Xu,u,Yq,Xq);
v_spaced = interp2(Yv,Xv,v,Yq,Xq);
end

clf
if FreeMyosin==1
    subaxis(1,3,1)
end
hold on
hfig = pcolor(Xm,Ym,m);
set(hfig,'EdgeColor','none')
colorbar
caxis([0,2])
if ActinFluid==1
hfig2 = quiver(Xq,Yq,u_spaced,v_spaced,'w');
end
for k = 1:length(Fibers)
    XYfiber = Fibers{k}{1};
    Zfiber = Fibers{k}{4};
    UVfiber = Fibers{k}{3}; Ufiber = UVfiber(:,1); Vfiber = UVfiber(:,2);
    if isempty(Zfiber) ~= 1
        plot(Zfiber(:,1),Zfiber(:,2),'r*','MarkerSize',3)
    end
    if Fibers{k}{8} == 1
    plot([XYfiber(:,1);XYfiber(1,1)],[XYfiber(:,2);XYfiber(1,2)],'w.-','LineWidth',3)
    else
    plot(XYfiber(:,1),XYfiber(:,2),'w.-','LineWidth',3)
    end
    quiver(XYfiber(1:1:end,1),XYfiber(1:1:end,2),Ufiber(1:1:end),Vfiber(1:1:end),'r')
end
ylabel('y','Fontsize',14)
xlabel('x','Fontsize',14)
axis([La Lb Lc Ld])
axis equal
if FreeMyosin==1
subaxis(1,3,2)
hold on
hfig = pcolor(Xm,Ym,b);
set(hfig,'EdgeColor','none')
colorbar
ylabel('y','Fontsize',14)
xlabel('x','Fontsize',14)
axis([La Lb Lc Ld])
axis equal

    subaxis(1,3,3)
    yyaxis('left')
    plot(m(:,round(size(m,2)/2)),'k')
    ylabel('Bound')
    yyaxis('right')
    plot(u(:,round(size(u,2)/2)),'r')
    % ylim([0,1])
    ylabel('velocity')
    title('Midplane myosin')


end

if ActinFluid==0
sgtitle(strcat('$t=',sprintf('%.2f',t),', \tilde{\mu}=$',sprintf('%.2f',eta2),',  $\zeta=$',sprintf('%.2f',zeta),',  $|U(x,y,t)|\leq$',sprintf('%.4f',max(sqrt(abs(Ufiber).^2+abs(Vfiber).^2),[],'all'))),'interpreter','latex','Fontsize',16)
else
sgtitle(strcat('$t=',sprintf('%.2f',t),', \tilde{\mu}=$',sprintf('%.2f',eta2),',  $\zeta=$',sprintf('%.2f',zeta),',  $|u(x,y,t)|\leq$',sprintf('%.4f',max(sqrt(abs(u_spaced).^2+abs(v_spaced).^2),[],'all'))),'interpreter','latex','Fontsize',16)
end

 pause(.01)
end
end

% stopping the simulation for any other reason
if stopSIM == 1
    disp('Simulation stopped')
    break
end


end


waitbar(1,wbar,'Finished!');
%pause(1)

close(wbar)

end
    


