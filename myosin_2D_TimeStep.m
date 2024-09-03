function [m] = myosin_2D_TimeStep(u,v,m,hx,hy,dt,source,AdvFlux_choice,BC_choice)

if BC_choice == 0
    u = [u;u(1,:)];
    v = [v,v(:,1)];
end

% calculate DCU flux at faces in interior
f_DCU = zeros(size(u));
f_DCU(2:end-1,:) = (u(2:end-1,:)>0).*u(2:end-1,:).*m(1:end-1,:) + (u(2:end-1,:)<0).*u(2:end-1,:).*m(2:end,:);
% flux at left boundary and right boundary % for now, periodic B.C.
if BC_choice==0
f_DCU(1,:) = (u(1,:)>0).*u(1,:).*m(end,:) + (u(1,:)<0).*u(1,:).*m(1,:);
f_DCU(end,:) = (u(end,:)>0).*u(end,:).*m(end,:) + (u(end,:)<0).*u(end,:).*m(1,:);
elseif BC_choice==1
f_DCU(1,:) = (u(1,:)<0).*u(1,:).*m(1,:);
f_DCU(end,:) = (u(end,:)>0).*u(end,:).*m(end,:);
end

g_DCU = zeros(size(v));
g_DCU(:,2:end-1) = (v(:,2:end-1)>0).*v(:,2:end-1).*m(:,1:end-1) + (v(:,2:end-1)<0).*v(:,2:end-1).*m(:,2:end);
% flux at top and bottom boundaries % for now, periodic B.C.
if BC_choice==0
g_DCU(:,1) = (v(:,1)>0).*v(:,1).*m(:,end) + (v(:,1)<0).*v(:,1).*m(:,1);
g_DCU(:,end) = (v(:,end)>0).*v(:,end).*m(:,end) + (v(:,end)<0).*v(:,end).*m(:,1);
elseif BC_choice==1
g_DCU(:,1) = (v(:,1)<0).*v(:,1).*m(:,1);
g_DCU(:,end) = (v(:,end)>0).*v(:,end).*m(:,end);
end

ff = f_DCU; gg = g_DCU;

if AdvFlux_choice==1
% calculate CTU fluxes at faces in interior
f_CTU = zeros(size(u));

Bpm =  zeros(size(v(:,1:end-1)));
Bpm(:,2:end) = (v(:,2:end-1)>0).*v(:,2:end-1).*(m(:,2:end)-m(:,1:end-1));
% at bottom, for now periodic B.C.
if BC_choice==0
Bpm(:,1) = (v(:,1)>0).*v(:,1).*(m(:,1)-m(:,end));
elseif BC_choice==1
Bpm(:,1) = (v(:,1)>0).*v(:,1).*(m(:,1));
end

Bmm = zeros(size(v(:,1:end-1)));
Bmm(:,1:end-1) = (v(:,2:end-1)<0).*v(:,2:end-1).*(m(:,2:end)-m(:,1:end-1));
% at top, for now periodic B.C.
if BC_choice==0
Bmm(:,end) = (v(:,end)<0).*v(:,end).*(m(:,1)-m(:,end));
elseif BC_choice==1
Bmm(:,end) = (v(:,end)<0).*v(:,end).*(-m(:,end));
end

f_CTU(2:end-1,:) = (u(2:end-1,:)>0).*u(2:end-1,:).*(Bpm(1:end-1,:) + Bmm(1:end-1,:)) + ...
                            (u(2:end-1,:)<0).*u(2:end-1,:).*(Bpm(2:end,:) + Bmm(2:end,:));
% periodic B.C.
if BC_choice==0
% left
f_CTU(1,:) = (u(1,:)>0).*u(1,:).*(Bpm(end,:)+Bmm(end,:)) + ...
                            (u(1,:)<0).*u(1,:).*(Bpm(1,:)+Bmm(1,:));
% right
f_CTU(end,:) = (u(end,:)>0).*u(end,:).*(Bpm(end,:)+Bmm(end,:)) + ...
                            (u(end,:)<0).*u(end,:).*(Bpm(1,:)+Bmm(1,:));
elseif BC_choice==1
% left
f_CTU(1,:) = (u(1,:)<0).*u(1,:).*(Bpm(1,:)+Bmm(1,:));
% right
f_CTU(end,:) = (u(end,:)>0).*u(end,:).*(Bpm(end,:)+Bmm(end,:));
end

f_CTU = -1/2*dt/hy*f_CTU;

g_CTU = zeros(size(v));

Apm = zeros(size(u(1:end-1,:)));
Apm(2:end,:) = (u(2:end-1,:)>0).*u(2:end-1,:).*(m(2:end,:)-m(1:end-1,:));
% left BC periodic
if BC_choice==0
Apm(1,:) = (u(1,:)>0).*u(1,:).*(m(1,:)-m(end,:));
elseif BC_choice==1
Apm(1,:) = (u(1,:)>0).*u(1,:).*(m(1,:));
end

Amm = zeros(size(u(1:end-1,:)));
Amm(1:end-1,:) = (u(2:end-1,:)<0).*u(2:end-1,:).*(m(2:end,:)-m(1:end-1,:));
% right BC periodic
if BC_choice==0
Amm(end,:) = (u(end,:)<0).*u(end,:).*(m(1,:)-m(end,:));
elseif BC_choice==1
Amm(end,:) = (u(end,:)<0).*u(end,:).*(-m(end,:));
end

g_CTU(:,2:end-1) = (v(:,2:end-1)>0).*v(:,2:end-1).*(Apm(:,1:end-1)+Amm(:,1:end-1)) + ...
                    (v(:,2:end-1)<0).*v(:,2:end-1).*(Apm(:,2:end)+Amm(:,2:end));
if  BC_choice==0
% bottom % periodic BC
g_CTU(:,1) = (v(:,1)>0).*v(:,1).*(Apm(:,end)+Amm(:,end)) + (v(:,1)<0).*v(:,1).*(Apm(:,1)+Amm(:,1));
% top
% g_CTU(:,end) = (v(:,end)>0).*v(:,end).*(Apm(:,end)+Amm(:,end)) + ...
%               (v(:,end)<0).*v(:,end).*(Apm(:,end)+Amm(:,end)); %old line with presummed bug
g_CTU(:,end) = (v(:,end)>0).*v(:,end).*(Apm(:,end)+Amm(:,end)) + (v(:,end)<0).*v(:,end).*(Apm(:,1)+Amm(:,1));
elseif BC_choice==1
% bottom
g_CTU(:,1) =  (v(:,1)<0).*v(:,1).*(Apm(:,1)+Amm(:,1));
% top
g_CTU(:,end) = (v(:,end)>0).*v(:,end).*(Apm(:,end)+Amm(:,end));
end


g_CTU = -1/2*dt/hx*g_CTU;

ff = ff + f_CTU; 
gg = gg + g_CTU; 
end

div_ff = 1/hx*(ff(2:end,:)-ff(1:end-1,:));
div_gg = 1/hy*(gg(:,2:end)-gg(:,1:end-1));

m = m - dt*div_ff - dt*div_gg + dt*source;


end