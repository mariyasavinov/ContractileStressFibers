function [Fiber1] = construct_FiberElement(XYfiber, hq_array, Nb, kfiber, kBfiber, kadh, xi, contrConst, FiberAdhesion,FiberBending,FiberViscous,FiberContractile,FiberLoop)

if FiberAdhesion==1      % adhesive points at beginning and end of the fiber
    Zfiber = [XYfiber(1,:),1; XYfiber(end,:),Nb];
elseif FiberAdhesion==2  % ahdesions along length of the fiber
    kadh = kadh*mean(hq_array); % scale adhesive force by number of points
    Zfiber = [XYfiber, linspace(1,length(XYfiber),length(XYfiber))'];
else
    Zfiber = [];
end

% UVfiber
UVfiber = [];

% matrices for timestepping
% B for bending forces
if FiberBending==1
[B] = build_B(Nb,hq_array(1),kBfiber,FiberLoop);
if FiberViscous==1 || FiberContractile==3 %%% if you are doing FPI
    % reorder so it is X1;Y1;X2;Y2;...
    B_reorder = sparse(2*Nb,2*Nb);
    B_reorder(1:2:end,1:2:end) = B(1:Nb,1:Nb);
    % B_reorder(2:2:end,1:2:end) = B(Nb+1:end,Nb+1:end);
    B_reorder(2:2:end,2:2:end) = B(Nb+1:end,Nb+1:end);
    B = B_reorder;
end
else
B = sparse(2*Nb,2*Nb);
end
% identity matrix
I2Nb = speye(2*Nb);

% C and D matrices for Slip with FPI
C = sparse(2*Nb,2*Nb); 
D = sparse(2*Nb,2*Nb);

% bundle fiber parameters
fiberparams = {kfiber, kBfiber, kadh, xi, contrConst};

% build Fiber element
Fiber1 = {XYfiber, hq_array, UVfiber, Zfiber, {B,C,D}, I2Nb, fiberparams, FiberLoop};