function F = viscousForce(XYfiber,Ufiber,Vfiber,hq_array,xi,FiberLoop)

% calculate viscous force due to fiber's configuration (along fiber)
Nb = size(XYfiber,1);

% tensions along fiber
% T = kfiber*(vecnorm(XYfiber(2:end,:)-XYfiber(1:end-1,:),2,2)./hq_array(1:Nb-1)'-1);

% unit vectors
tau = (XYfiber(2:end,:)-XYfiber(1:end-1,:))./vecnorm(XYfiber(2:end,:)-XYfiber(1:end-1,:),2,2);
Tu = xi*((Ufiber(2:end)-Ufiber(1:end-1)).*tau(:,1)+(Vfiber(2:end)-Vfiber(1:end-1)).*tau(:,2))./hq_array(1:Nb-1)';

B = Tu.*tau; B = [B;[0,0]];
% if the fiber is connected
if FiberLoop == 1
tau_conn = (XYfiber(1,:)-XYfiber(end,:))./vecnorm(XYfiber(1,:)-XYfiber(end,:),2,2);
Tu_conn = xi*((Ufiber(1)-Ufiber(end)).*tau_conn(1)+(Vfiber(1)-Vfiber(end)).*tau_conn(2))/hq_array(Nb);
B(end,:) = Tu_conn.*tau_conn;
end

% force on all fiber points
F = (B - [B(end,:);B(1:end-1,:)]);

end