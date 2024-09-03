function [Ufiber] = u_interp1D(Xu,Yu,XYfiber,u)

Nb = size(XYfiber,1);
hx = Xu(2,1)-Xu(1,1);
hy = Yu(1,2)-Yu(1,1);

% interpolate velocity field to velocity of fiber points
Ufiber = zeros(Nb,1);

% for U component, get body position relative to Xu,Yu grid
s = (XYfiber-[Xu(1,1),Yu(1,1)]).*[1/hx,1/hy];
i = floor(s);
r = s-i;

phi1 = phi_vec(r(:,1),1);
phi2 = phi_vec(r(:,2),2);

% evaluate delta function
w = phi1.*phi2;
w = permute(w,[1,3,2]);

for k = 1:Nb
   
    % find adjacent fluid points where u component of velocity is defined
    i1 = (i(k,1)):(i(k,1)+3); %mod((i(k,1)-1):(i(k,1)+2),Nx+2)+1; for periodic
    i2 = (i(k,2)):(i(k,2)+3); %mod((i(k,2)-1):(i(k,2)+2),Ny+1)+1; for perioidc 
   
    % interpolate fluid velocity in x direction
    ww = w(:,:,k);
    Uk = sum(sum(ww.*u(i1,i2)));
    Ufiber(k) = Uk;
   
end

end