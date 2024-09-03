function F = contractileForce(XYfiber,Ufiber,Vfiber,hq_array,contrConst,FiberLoop,FiberContractile)

% calculate elastic force due to fiber's configuration (along fiber)

Nb = size(XYfiber,1);

l0sarc = contrConst(1); Fsarc = contrConst(2); Fs = contrConst(3); v0 = contrConst(4);
Nb = size(XYfiber,1);
tau = (XYfiber(2:end,:)-XYfiber(1:end-1,:))./vecnorm(XYfiber(2:end,:)-XYfiber(1:end-1,:),2,2);
if FiberContractile==1              % = F*N_s = F*hq/l0sarc;
    Tcontr = Fsarc*hq_array(1:Nb-1)'/l0sarc;
elseif FiberContractile==2          % = F*N_s = F*l/l0sarc;
    Tcontr = Fsarc*vecnorm(XYfiber(2:end,:)-XYfiber(1:end-1,:),2,2)/l0sarc;
elseif FiberContractile==3          % = Fs*(1-v/v0)
    Tcontr = Fs*(1+1/v0*(((Ufiber(2:end)-Ufiber(1:end-1)).*tau(:,1)+(Vfiber(2:end)-Vfiber(1:end-1)).*tau(:,2)))./hq_array(1:Nb-1)');
end
B = Tcontr.*tau; B = [B;[0,0]];
if FiberLoop == 1
tau_conn = (XYfiber(1,:)-XYfiber(end,:))./vecnorm(XYfiber(1,:)-XYfiber(end,:),2,2);
if FiberContractile==1              % = F*N_s = F*hq/l0sarc;
    Tcontr_conn = Fsarc*hq_array(Nb)/l0sarc;
elseif FiberContractile==2          % = F*N_s = F*l/l0sarc;
    Tcontr_conn = Fsarc*vecnorm(XYfiber(1,:)-XYfiber(end,:),2,2)/l0sarc;
elseif FiberContractile==3          % = Fs*(1-v/v0)
    Tcontr_conn = Fs*(1+1/(v0)*(((Ufiber(1)-Ufiber(end)).*tau_conn(1)+(Vfiber(1)-Vfiber(end)).*tau_conn(2)))/hq_array(Nb));
end
B(end,:) = Tcontr_conn.*tau_conn;
end
% force on all fiber points
F = (B - [B(end,:);B(1:end-1,:)]);

end