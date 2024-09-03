function [Ffiber_spread, felastic, fbending, fviscous, fadhesion, fcontractile] = get_Ffiber_spread(XYuv,XYfiber,Ufiber,Vfiber,Zfiber,hq_array,kfiber,kBfiber,xi,contrConst,kadh,BC_choice,Fiber_vars)

FiberLoop = Fiber_vars(1); FiberContractile = Fiber_vars(2); FiberAdhesion = Fiber_vars(3);
FiberViscous = Fiber_vars(4); FiberBending = Fiber_vars(6);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elastic force on fiber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
felastic = elasticForce(XYfiber,hq_array,kfiber,FiberLoop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bending force on fiber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fbending = zeros(size(felastic));
if FiberBending==1
    fbending = bendingForce(XYfiber,hq_array,kBfiber,FiberLoop);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% viscous force on fiber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fviscous = zeros(size(felastic));
if FiberViscous==1
    fviscous = viscousForce(XYfiber,Ufiber,Vfiber,hq_array,xi,FiberLoop);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adhesion forces on fiber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fadhesion = zeros(size(felastic));
% if FiberAdhesion == 1
if isempty(Zfiber)~=1
    for k = 1:size(Zfiber,1)
        fadhesion(Zfiber(k,3),:) = fadhesion(Zfiber(k,3),:)-kadh*(XYfiber(Zfiber(k,3),:)-Zfiber(k,1:2));
        % if Zfiber(k,3)==1 || Zfiber(k,3)==length(XYfiber)  %% TEMPORARY
        %     fadhesion(Zfiber(k,3),:) = fadhesion(Zfiber(k,3),:)-kadh*(XYfiber(Zfiber(k,3),:)-Zfiber(k,1:2));
        % else
        %     fadhesion(Zfiber(k,3),:) = fadhesion(Zfiber(k,3),:)-kadh/4*(XYfiber(Zfiber(k,3),:)-Zfiber(k,1:2));
        % end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contractile forces on fiber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcontractile = zeros(size(felastic));
% l0sarc = contrConst(1); Fsarc = contrConst(2); Fs = contrConst(3); v0 = contrConst(4);
if FiberContractile~=0
    fcontractile = contractileForce(XYfiber,Ufiber,Vfiber,hq_array,contrConst,FiberLoop,FiberContractile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total sum of fiber forces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ffiber = felastic + fbending + fviscous + fadhesion + fcontractile;

% spread F to the Xu,Yu grid (for x-component of f)
fu = fu_spread1D(Xu,Yu,XYfiber,Ffiber);
fv = fu_spread1D(Xv,Yv,XYfiber,Ffiber);
fu = fu(:,:,1); % only x-direction of force for u
fv = fv(:,:,2); % only y-direction of force for v
% reshape into vector
if BC_choice==0
    fu = reshape(fu,Nx*Ny,1);
    fv = reshape(fv,Nx*Ny,1);
elseif BC_choice==1
    fu = reshape(fu,(Nx+2)*(Ny+1),1);
    fv = reshape(fv,(Nx+1)*(Ny+2),1);
elseif BC_choice==2
    fu = reshape(fu,(Nx+1)*(Ny+1),1);
    fv = reshape(fv,(Nx+1)*(Ny+2),1);
end
Ffiber_spread = [fu; fv];

end