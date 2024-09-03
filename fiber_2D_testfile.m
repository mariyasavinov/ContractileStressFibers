
% author: Mariya Savinov

clear all
close all
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests / Plotting choices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PlotWhileRun = 1;               % = 1 if you want solver to plot every 0.05 time units

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FixedPointIteration = 0; % unless other quantities set it differently
FPItol = 1e-5; %1e-6

ActinFluid = 1; % = 0 if no actin fluid
                % = 1 if with an actin fluid

BC_choice = 1; % = 0 if periodic
               % = 1 if no-stress
               % = 2 if periodic in x
               
IC_myosin = 1; % = 0 if no myosin
               % = 1 if myosin is a gaussian peak at center of domain
               % = 2 if myosin is a gaussian peak off center diagonally
               % = 3 if myosin is two gaussian peaks, both off-center
               % = 4 if myosin is constant = 1 over entire domain
               % = 5 if myosin is two gaussian peaks, both on right side in corners
               % = 6 f myosin is (sine(x)*sin(y))^6 peak at center of domain

IC_b = 3;      % as above for myosin

IC_fiber = 41;  % = 0 if circle with radius r_deformed at center of domain with r_rest rest configuration size
               % = 1 if a vertical line with length L_deformed at center of domain with L_rest rest configuration length 
               % = 2 if random orientation at rest with length 0.5
               % = 3 if stretched to curve from rest length = 0.5
               % = 41,42,43,44 for IC near edges of the domain
               % = 5 for horizontal line with length L_deformed at center of domain with L_rest rest configuration length
               % = 6 if stretched to curve from rest length = 0.5, which is then set to be the rest configuration 
               
FiberLoop = 0; % = 1 if fiber is a connected loop

FiberBending = 0; % = 1 if fiber has forces from bending

FiberContractile = 1;
               % = 0 if no contraction
               % = 1 if fixed contractile force Fm = F*N_s, where N_s = number of unstreched sarcomere units 
               % = 2 if length-dependent contractile force Fm = F*N_s, where N_s = number of unstretched sarcomere units that could fit in current length
               % = 3 if velocity-dependent contractile force Fm = Fs(1-1/v0*dl/dt), where Fs is the stall force and v0 is the free velocity 

               if FiberContractile==3
                   FixedPointIteration = 1;
               end

FiberViscous = 1;
               % = 1 if fiber is viscoelastic

               if FiberViscous==1
                   FixedPointIteration = 1;
               end

FiberAdhesion = 0;
               % = 1 if fiber is adhering to target points at ends
               % = 2 if fiber is adhesive all along length

FiberSlip = 0;
               % = 0 if no-slip BC between fiber and fluid, = 1 if slip

if FiberSlip == 1 && ActinFluid == 1
    error('For 2D case, cannot do slip and IBM together. Choose either a drag model case or IBM case.')
end

if FiberSlip == 0 && ActinFluid == 0
    error('Fiber will not move. Set either FiberSlip==1 or ActinFluid==1')
end

LaserAblation = 0;
               % = 1 if laser ablation at Indx_Laser at time T_Laser
               
BC_order = 2;  % = 1 if first order
               % = 21 if second order except at corners
               % = 2 if second order everywhere
             
if (IC_fiber == 1 || IC_fiber == 2) && FiberLoop == 1
    error('To generate a connected loop fiber, pick a closed rest configuration')
end


if FiberSlip==1 && FiberViscous==1 && FiberLoop==1
    error('Viscous force matrix not implemented for loop fiber')
end

AdvFlux_choice = 1;
            % = 0 for DCU flux
            % = 1 for CTU flux

FreeMyosin = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st viscosity coefficient
eta1 = 0.6;                      % nondimensional = (hydrodynamic length / domain lengthscale)^2
% 2nd viscosity coefficient
eta2 = 0.6;                      % must at minimum be >= 1/3*eta1 by thermodynamics arguments
% drag coefficient
zeta = 1.0;                      % Do not change from =1 (due to non-dimensionalization)
% contraction coefficient 
k = 1.0;                         % Do not change from =1 (due to non-dimensionalization)

% myosin exchange bound --> free
alpha = 1.0;
% myosin exchange free --> bound
beta = 1.0;
% free myosin diffusion constant
d0 = 0.001;


% fiber force scale
f0 = 1.0;
% fiber elastic strength/modulus
kfiber = 1;
% fiber bending strength/modulus
kBfiber = 0.0005;
% fiber adhesion strengths (note: kadh = k_adh*deltaZ, so measure of adhesion is absorbed into this constant)
kadh = 200;

% fiber viscosity coefficient
xi = 0.05;  %0.25 may be maximum without slip term

% fiber contractility constants
l0sarc = 0.024;             % rescaled unstretched sarcomere length, assuming to be ~ 2.4 microns
Fsarc = 0.5;                % contractile force from a single sarcomere
Fs = 0.025;                   % stall force
v0 = 15;                   % free velocity per unit length
contrConst = [l0sarc,Fsarc,Fs,v0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain and Spatial Discretization Sizes
% e.g., consider stress fiber is 0.5 microns, domain is R = 100 microns. h~0.42 microns, R_h = 1.2*h = 0.5 microns 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain [La,Lb]x[Lc,Ld]
La = 0; Lb = 1; Lc = 0; Ld = 1;
% fluid discretization sizes
Nx=100; Ny = 100;
% Nx = 256; Ny = 256;
% Nx = 240; Ny = 240;
% fiber discretization size
Nb = 60;

if ActinFluid==1
    % fiber slip boundary condition coefficient 
    Rh = 1.2/Nx;        % approximate hydrodynamic radius from IBM
    Rsf = 0.5/100;      % approximate ratio of stress fiber thickness to domain size 0.5 microns / 100 microns 
    gamma = Rh*Rsf/(Rh-Rsf)*12*pi*eta1*(eta2+eta1)*(3*eta2-eta1)/(eta2+3*eta1)/(2*eta2+3*eta1);
else
    gamma = 1.0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% laser ablation protocol1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Indx_Laser = round(Nb/2);          % break between Indx_Laser and Indx_Laser+1
T_Laser = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timestepping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfinal = 15; 
dt = 0.001; % 0.0001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save vars to pass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global_vars = [ActinFluid, BC_choice, BC_order, PlotWhileRun, FixedPointIteration, FPItol, AdvFlux_choice,FreeMyosin];
Fiber_vars = [FiberLoop, FiberContractile, FiberAdhesion, FiberViscous, FiberSlip, FiberBending, LaserAblation, Indx_Laser, T_Laser];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resting and initial deformed fiber lengths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_rest = 0.25;      % only for circle
r_deformed = 0.25;  % only for circle
L_rest = 0.5; 
L_deformed = 0.5;%0.75;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bound and free myosin initial distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m_func, ~, ~] = myosin_ICfunc(IC_myosin,La,Lb,Lc,Ld);
[b_func, ~, ~] = myosin_ICfunc(IC_b,La,Lb,Lc,Ld);

%% Solution test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build the Eulerian grid, get grid points for u and v, get grid sizes
[X,Y,Xu,Yu,Xv,Yv,hx,hy] = build_EulerianGrid(Nx,Ny,La,Lb,Lc,Ld,BC_choice);
if BC_choice ==0
    Xm = X(1:end,1:end)+hx/2;
    Ym = Y(1:end,1:end)+hy/2;
else
    Xm = X(1:end-1,1:end-1)+hx/2;
    Ym = Y(1:end-1,1:end-1)+hy/2;
end
% build cell with all grid points to pass to solver
XYuvm = {X,Y,Xu,Yu,Xv,Yv,Xm,Ym};

% initial myosin configuration
m0 = m_func(Xm,Ym);
b0 = b_func(Xm,Ym);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build fibers (Lagrangian)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IC_fiber~=45 %1 && IC_fiber~=42 && IC_fiber~=43 && IC_fiber~=44
% first fiber
[XYfiber_rest, hq_array, XYfiber] = fiber_IC(IC_fiber,r_rest,r_deformed,L_rest,L_deformed,Nb,La,Lb,Lc,Ld,FiberLoop);
Nb = size(XYfiber,1);
Fiber1 = construct_FiberElement(XYfiber, hq_array, Nb, kfiber, kBfiber, kadh, xi, contrConst, FiberAdhesion,FiberBending,FiberViscous,FiberContractile,FiberLoop);
Fibers = {Fiber1};

else  % 4 fibers at edges of domain 

% first fiber
[XYfiber_rest, hq_array, XYfiber] = fiber_IC(41,r_rest,r_deformed,L_rest,L_deformed,Nb,La,Lb,Lc,Ld,FiberLoop);
Nb = size(XYfiber,1);
Fiber1 = construct_FiberElement(XYfiber, hq_array, Nb, kfiber, kBfiber, kadh, xi, contrConst, FiberAdhesion,FiberBending,FiberViscous,FiberContractile,FiberLoop);

% second fiber
[XYfiber_rest, hq_array, XYfiber] = fiber_IC(42,r_rest,r_deformed,L_rest,L_deformed,Nb,La,Lb,Lc,Ld,FiberLoop);
Nb = size(XYfiber,1);
Fiber2 = construct_FiberElement(XYfiber, hq_array, Nb, kfiber, kBfiber, kadh, xi, contrConst, FiberAdhesion,FiberBending,FiberViscous,FiberContractile,FiberLoop);

% third fiber
[XYfiber_rest, hq_array, XYfiber] = fiber_IC(43,r_rest,r_deformed,L_rest,L_deformed,Nb,La,Lb,Lc,Ld,FiberLoop);
Nb = size(XYfiber,1);
Fiber3 = construct_FiberElement(XYfiber, hq_array, Nb, kfiber, kBfiber, kadh, xi, contrConst, FiberAdhesion,FiberBending,FiberViscous,FiberContractile,FiberLoop);

% fourth fiber
[XYfiber_rest, hq_array, XYfiber] = fiber_IC(44,r_rest,r_deformed,L_rest,L_deformed,Nb,La,Lb,Lc,Ld,FiberLoop);
Nb = size(XYfiber,1);
Fiber4 = construct_FiberElement(XYfiber, hq_array, Nb, kfiber, kBfiber, kadh, xi, contrConst, FiberAdhesion,FiberBending,FiberViscous,FiberContractile,FiberLoop);

Fibers = {Fiber3, Fiber1, Fiber2, Fiber4};
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot initial configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FreeMyosin~=1
figIC = figure('position',[50,50,500,500]);
else
figIC = figure('position',[50,50,1100,500]);
subaxis(1,2,1)
end
hfig2 = pcolor(Xm,Ym,m0);
hold on
set(hfig2,'EdgeColor','none')
colorbar
for FiberIndx = 1:length(Fibers)
    XYfiber = Fibers{FiberIndx}{1};
    Zfiber = Fibers{FiberIndx}{4};
    if Fibers{FiberIndx}{8} == 1
    plot([XYfiber(:,1);XYfiber(1,1)],[XYfiber(:,2);XYfiber(1,2)],'w.-','LineWidth',3)
    else
    plot(XYfiber(:,1),XYfiber(:,2),'w.-','LineWidth',3)
    end
    if isempty(Zfiber) ~= 1
        plot(Zfiber(:,1),Zfiber(:,2),'r*','MarkerSize',3)
    end
end
ylabel('y','Fontsize',20)
xlabel('x','Fontsize',20)
axis([La Lb Lc Ld])
axis equal


if FreeMyosin==1
title('bound myosin','FontSize',13)
subaxis(1,2,2)
hfig = pcolor(Xm,Ym,b0);
set(hfig,'EdgeColor','none')
colorbar
ylabel('y','Fontsize',14)
xlabel('x','Fontsize',14)
clim([0,1])
axis([La Lb Lc Ld])
axis equal
title('free myosin','FontSize',13)
end

sgtitle('Initial Condition')
pause(1)
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsteps = round(tfinal/dt);
% solve for fluid and fiber velocities + fiber position for nsteps
[Fibers_array,u_array,v_array,m_array,b_array,t_array,iterK0] =  fiber_2D_solver(XYuvm,hx,hy,Fibers,nsteps,dt,eta1,eta2,zeta,d0,alpha,beta,gamma,f0,m0,b0,global_vars,Fiber_vars);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct Extra Arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%
% get array length over time
L0_array = [];
L_array = [];
for tindx=1:length(t_array)
    for Fiberk=1 %k=1:length(Fibers_XYfiber_array{1})
        XYfiber = Fibers_array{tindx}{Fiberk}{1};
        L0_array(tindx) = sum(Fibers_array{tindx}{Fiberk}{2});
        % get rest and current length of fiber
        L = 0;
        for j=1:length(XYfiber)-1
            L = L + norm(XYfiber(j+1,:)-XYfiber(j,:));
        end
        L_array(tindx) = L;
    end
end

% get curvature over time
Cmean_array = [];
Cmax_array = [];
Fiberk=1;
for tindx = 1:length(t_array)
    XYfiber = Fibers_array{tindx}{Fiberk}{1};
    hq = Fibers_array{tindx}{Fiberk}{2}(1);
    K = (XYfiber(3:Nb,:) - 2*XYfiber(2:Nb-1,:) + XYfiber(1:Nb-2,:))/hq^2;
    K = vecnorm(K,2,2);
    Cmean_array(tindx) = mean(K);
    Cmax_array(tindx) = max(K);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%
% Plot initial and final configuration
arrow_num = 100;
x_spaced = La:(Lb-La)/arrow_num:Lb;
y_spaced = Lc:(Ld-Lc)/arrow_num:Ld;
[Xq,Yq] = meshgrid(x_spaced,y_spaced);

if ActinFluid==1
u_max = 0;
U_max = 0;
for m=1:length(u_array)-1
    u_spaced = interp2(Yu,Xu,u_array{m},Yq,Xq);
    v_spaced = interp2(Yv,Xv,v_array{m},Yq,Xq);
    u_norm = max(sqrt(u_spaced.^2+v_spaced.^2),[],'all');
    u_max = max(u_max,u_norm);
end
arrow_num = 30;
x_spaced = La:(Lb-La)/arrow_num:Lb;
y_spaced = Lc:(Ld-Lc)/arrow_num:Ld;
[Xq,Yq] = meshgrid(x_spaced,y_spaced);
quivscale = hx/u_max;
quivscale = hx/u_max*20;
else
U_max = 0;
for tindx=1:length(t_array)-1
    UVfiber = Fibers_array{tindx}{1}{3}; Ufiber = UVfiber(:,1); Vfiber = UVfiber(:,2);
    U_norm = max(sqrt(UVfiber.^2+Vfiber.^2),[],'all');
    U_max = max(U_max,U_norm);
end
[Xq,Yq] = meshgrid(x_spaced,y_spaced);
quivscale = hx/U_max;
quivscale = hx/U_max*20;
end

for tindx = [1, round(length(t_array)/4),round(length(t_array)/2), round(3*length(t_array)/4), length(t_array)]
m = m_array{tindx};
if ActinFluid==1
u_spaced = interp2(Yu,Xu,u_array{tindx},Yq,Xq); v_spaced = interp2(Yv,Xv,v_array{tindx},Yq,Xq);
end
t = t_array(tindx);

fig1 = figure('position',[50,50,600,600]);
set(gca,'FontSize',14)
hold on
hfig2 = pcolor(Xm,Ym,m);
set(hfig2,'EdgeColor','none')
caxis([0,1])
colorbar
if ActinFluid==1; hfig2 = quiver(Xq,Yq,u_spaced*quivscale,v_spaced*quivscale,'w','AutoScale','off'); end

for k = 1:length(Fibers_array{tindx})
    XYfiber = Fibers_array{tindx}{k}{1};
    UVfiber = Fibers_array{tindx}{k}{3}; Ufiber = UVfiber(:,1); Vfiber = UVfiber(:,2);
    plot([XYfiber(:,1)],[XYfiber(:,2)],'r.-','LineWidth',3)
    % quiver(XYfiber(2:2:end,1),XYfiber(2:2:end,2),Ufiber(2:2:end)*quivscale,Vfiber(2:2:end)*quivscale,'r','AutoScale','off')
    quiver(XYfiber(2:2:end,1),XYfiber(2:2:end,2),Ufiber(2:2:end),Vfiber(2:2:end),'r')
end

axis([La Lb Lc Ld])
axis equal
title(strcat('$t=',sprintf('%.2f',t),', \tilde{\mu}=$',sprintf('%.2f',eta2),',  $\zeta=$',sprintf('%.2f',zeta)),'interpreter','latex','Fontsize',17)

end

%% Plot length L of fiber over time
% figure;
fig20 = figure('position',[50,50,500,500]);
set(gca,'FontSize',14)
hold on
plot(t_array,L0_array,'k.-','DisplayName','rest')
plot(t_array,L_array,'r.-','LineWidth',3,'DisplayName','L(t)')


%% Plot curvature measure over time
% figure;
fig50 = figure('position',[50,50,500,500]);
set(gca,'FontSize',14)
hold on
plot(t_array,Cmean_array,'k.-','DisplayName','Mean curvature')
plot(t_array,Cmax_array,'r.-','LineWidth',3,'DisplayName','Max curvature')

xlabel('time','Fontsize',16)
ylabel('|C(q)|','Fontsize',16)

g = fittype('b*exp(-c*x)');
func0 = fit(t_array(1:end)',Cmax_array(1:end)',g);
t_array_fine = linspace(t_array(1),t_array(end),2000);
fig50;
plot(t_array_fine,func0(t_array_fine),'b--','LineWidth',2,'DisplayName',sprintf('%.4f exp(-%.4f*x)',func0.b,func0.c));%sprintf('%.4f+%.4f exp(-%.4f*x)',func0.a,func0.b,func0.c));

title('Curvature of fiber','Fontsize',14)
legend('Fontsize',12,'location','best')




%% Plot velocity profile along fiber over time

fig3=figure('position',[50,50,500,500]);
set(gca,'FontSize',14)
hold on

for tindx = [5, round(length(t_array)/5), round(length(t_array)/5*2),...
        round(length(t_array)/5*3),round(length(t_array)/5*4),length(t_array)]

XYfiber = Fibers_array{tindx}{1}{1};
tau = (XYfiber(2:end,:)-XYfiber(1:end-1,:))./vecnorm(XYfiber(2:end,:)-XYfiber(1:end-1,:),2,2);
UVfiber = Fibers_array{tindx}{1}{3}; Ufiber = UVfiber(:,1); Vfiber = UVfiber(:,2);

velocity = Ufiber(1:end-1).*tau(:,1)+Vfiber(1:end-1).*tau(:,2);
t = t_array(tindx);

plot(velocity,'o-','LineWidth',3,'DisplayName',sprintf('t=%.2f',t))

end
xlabel('index along fiber')
ylabel('fiber velocity')
legend('location','best','FontSize',12)

