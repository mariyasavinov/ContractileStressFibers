function [XYfiber_rest, hq_array, XYfiber] = fiber_IC(IC_fiber,r_rest,r_deformed,L_rest,L_deformed,Nb,La,Lb,Lc,Ld,FiberLoop)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resting and initial deformed fiber configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IC_fiber == 0                 % circle
center_point = [(Lb+La)/2,(Ld+Lc)/2];
elseif IC_fiber == 1             % vertical line
center_point = [(Lb+La)*1/2,(Ld+Lc)*1/2];
elseif IC_fiber == 2             % random
starting_point = [La,Lc]+[round(rand(1)*rand(1),4),round(rand(1)*rand(1),4)];
initial_angle = rand(1)*pi/2;
elseif IC_fiber == 3
starting_point = [(Lb+La)/2-L_rest/2,(Ld+Lc)/2];
% end_point = [(Lb+La)/2+L_rest/2,(Ld+Lc)/2];
elseif IC_fiber == 41  
starting_point = [La+0.125*(Lb-La),Lc+0.1*(Ld-Lc)];
end_point = [Lb-0.125*(Lb-La),Lc+0.1*(Ld-Lc)];
elseif IC_fiber == 42
starting_point = [La+0.1*(Lb-La),Lc+0.125*(Ld-Lc)];
end_point = [La+0.1*(Lb-La),Ld-0.125*(Ld-Lc)];
elseif IC_fiber == 43  
starting_point = [La+0.125*(Lb-La),Ld-0.1*(Ld-Lc)];
end_point = [Lb-0.125*(Lb-La),Ld-0.1*(Ld-Lc)];
elseif IC_fiber == 44  
starting_point = [Lb-0.1*(Lb-La),Lc+0.125*(Ld-Lc)];
end_point = [Lb-0.1*(Lb-La),Ld-0.125*(Ld-Lc)];
elseif IC_fiber == 5
center_point = [(Lb+La)*1/2,(Ld+Lc)*1/2];
elseif IC_fiber == 6
starting_point = [(Lb+La)/2-L_rest/2,(Ld+Lc)/2];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build rest configuration and get hq_array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IC_fiber == 0
% rest configuration positions
XYfiber_rest = [center_point(1)+r_rest*cos(2*pi/(Nb)*[0:Nb-1]);center_point(2)+r_rest*sin(2*pi/(Nb)*[0:Nb-1])]';
elseif IC_fiber == 1
% rest configuration positions
XYfiber_rest = [zeros(Nb,1)+center_point(1), (-L_rest/2:L_rest/(Nb-1):L_rest/2)'+center_point(2)];
elseif IC_fiber == 2
fiber_fits = 0;
mystep = 0.5/(Nb-1);
while fiber_fits == 0
    XYfiber_rest = zeros(Nb,2);
    XYfiber_rest(1,:) = [starting_point];
    myangle = initial_angle;
    for kcount=1:Nb-1
        XYfiber_rest(kcount+1,:) = XYfiber_rest(kcount,:)+mystep*[cos(myangle),sin(myangle)];
        myangle = myangle+0.1*(rand(1)-0.5)*pi/2;
    end

    if any(XYfiber_rest(:,1)<La) || any(XYfiber_rest(:,1)>Lb) || any(XYfiber_rest(:,2)<Lc) || any(XYfiber_rest(:,2)>Ld)
        fiber_fits=0;
        starting_point = [La,Lc]+[round(rand(1)*rand(1),4),round(rand(1)*rand(1),4)];
        initial_angle = rand(1)*pi/2;
    else
        fiber_fits=1;
    end
end
elseif IC_fiber == 3
    XYfiber_rest = [(0:L_rest/(Nb-1):L_rest)'+starting_point(1),zeros(Nb,1)+starting_point(2)];
elseif IC_fiber == 41 || IC_fiber == 42 || IC_fiber == 43 || IC_fiber == 44
    XYfiber_rest = [linspace(starting_point(1),end_point(1),Nb)',linspace(starting_point(2),end_point(2),Nb)'];
elseif IC_fiber == 5
% rest configuration positions
XYfiber_rest = [(-L_rest/2:L_rest/(Nb-1):L_rest/2)'+center_point(1), zeros(Nb,1)+center_point(2)];
elseif IC_fiber == 6

x0 = L_rest;
y0 = 0;
s0 = L_deformed;

a = y0/x0+0.1;

for k=1:100
    u = y0/x0 + a*x0;
    ell = y0/x0 - a*x0;
    dIu = sqrt(1+u^2);
    Iu = 0.5 * (u*dIu+ log(u+dIu));
    dIell = sqrt(1+ell^2);
    Iell = 0.5 * (ell*dIell+ log(ell+dIell));
    s = 0.5*(Iu-Iell)/a;
    ds = 0.5*(a*x0*(dIu+dIell)+Iell-Iu)/a^2;
    da = (s-s0)/ds;
    a = a-da;
    if abs(da)<1e-6
        break
    end
end
b = y0/x0 - a*x0;

x=linspace(0,L_rest,1000);
p = zeros(length(x),2);
p(:,1) = x+starting_point(1);
p(:,2) = a*x.^2+b*x+starting_point(2);
XYfiber_rest = curvspace(p,Nb);

% x = linspace(starting_point(1),0.5,Nb);
% y = starting_point(2)+1/4*sin(x/0.5*pi);
% % z = linspace(0,10,length(x));
% % N = 50;
% p = [x',y'];%,z'];
% XYfiber_rest = curvspace(p,Nb);
end

% get hq spring lengths
hq_array = vecnorm(XYfiber_rest(2:end,:)'-XYfiber_rest(1:end-1,:)');
if FiberLoop==1
   hq_array = [hq_array,norm(XYfiber_rest(1,:)-XYfiber_rest(end,:))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build starting configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IC_fiber == 0
XYfiber = (XYfiber_rest-[ones(Nb,1)*center_point(1), ones(Nb,1)*center_point(2)])*r_deformed/r_rest + [ones(Nb,1)*center_point(1), ones(Nb,1)*center_point(2)];
elseif IC_fiber == 1
XYfiber = (XYfiber_rest-[ones(Nb,1)*center_point(1), ones(Nb,1)*center_point(2)])*L_deformed/L_rest + [ones(Nb,1)*center_point(1), ones(Nb,1)*center_point(2)];
elseif IC_fiber == 2
XYfiber = XYfiber_rest;
elseif IC_fiber == 41 || IC_fiber == 42 || IC_fiber == 43 || IC_fiber == 44
XYfiber = XYfiber_rest;
elseif IC_fiber == 5
XYfiber = (XYfiber_rest-[ones(Nb,1)*center_point(1), ones(Nb,1)*center_point(2)])*L_deformed/L_rest + [ones(Nb,1)*center_point(1), ones(Nb,1)*center_point(2)];
elseif IC_fiber == 6
XYfiber = XYfiber_rest;

% myTempXY = [linspace(0.25,0.75,Nb)',zeros(Nb,1)];
% hq_array = vecnorm(myTempXY(2:end,:)'-myTempXY(1:end-1,:)');

elseif IC_fiber == 3
% finding quadratic passing through starting_point and endpoint
% y = a*x^2+b*x from (0,y0) to (L_rest,y0)  [will shift later]
% Using Newton Raphson Method

x0 = L_rest;
y0 = 0;
s0 = L_deformed;

a = y0/x0+0.1;

for k=1:100
    u = y0/x0 + a*x0;
    ell = y0/x0 - a*x0;
    dIu = sqrt(1+u^2);
    Iu = 0.5 * (u*dIu+ log(u+dIu));
    dIell = sqrt(1+ell^2);
    Iell = 0.5 * (ell*dIell+ log(ell+dIell));
    s = 0.5*(Iu-Iell)/a;
    ds = 0.5*(a*x0*(dIu+dIell)+Iell-Iu)/a^2;
    da = (s-s0)/ds;
    a = a-da;
    if abs(da)<1e-6
        break
    end
end
b = y0/x0 - a*x0;

x=linspace(0,L_rest,1000);
p = zeros(length(x),2);
p(:,1) = x+starting_point(1);
p(:,2) = a*x.^2+b*x+starting_point(2);
XYfiber = curvspace(p,Nb);

end

