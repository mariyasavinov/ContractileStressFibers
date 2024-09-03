function F = elasticForce(XYfiber,hq_array,kfiber,FiberLoop)

% calculate elastic force due to fiber's configuration (along fiber)

Nb = size(XYfiber,1);

% tensions along fiber
T = kfiber*(vecnorm(XYfiber(2:end,:)-XYfiber(1:end-1,:),2,2)./hq_array(1:Nb-1)'-1);
tau = (XYfiber(2:end,:)-XYfiber(1:end-1,:))./vecnorm(XYfiber(2:end,:)-XYfiber(1:end-1,:),2,2);


B = T.*tau; B = [B;[0,0]];
% if the fiber is connected
if FiberLoop == 1
T_conn = kfiber*(vecnorm(XYfiber(1,:)-XYfiber(end,:),2,2)/hq_array(Nb)-1);
tau_conn = (XYfiber(1,:)-XYfiber(end,:))./vecnorm(XYfiber(1,:)-XYfiber(end,:),2,2);
B(end,:) = T_conn.*tau_conn;
end

% force on all fiber points
F = (B - [B(end,:);B(1:end-1,:)]);

end