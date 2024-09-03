function fu = fu_spread1D(Xu,Yu,XYfiber,F)
% spread F to the Xu,Yu grid (for x-component of f)

Nb = size(XYfiber,1);
hx = Xu(2,1)-Xu(1,1);
hy = Yu(1,2)-Yu(1,1);
Nxu = size(Xu,1);
Nyu = size(Xu,2);

c = 1/(hx*hy);

fu = zeros(Nxu,Nyu,2);

% for fu component, get body position relative to Xu,Yu grid
s = (XYfiber-[Xu(1,1),Yu(1,1)]).*[1/hx,1/hy];
i = floor(s);
r = s-i;

phi1 = phi_vec(r(:,1),1);
phi2 = phi_vec(r(:,2),2);

% evaluate delta function
w = phi1.*phi2;
w = permute(w,[1,3,2]);

for k = 1:Nb
   
    % find adjacent fluid points which to which fiber force are speared
    i1 = (i(k,1)):(i(k,1)+3); %mod((i(k,1)-1):(i(k,1)+2),Nx+2)+1; for periodic
    i2 = (i(k,2)):(i(k,2)+3); %mod((i(k,2)-1):(i(k,2)+2),Ny+1)+1; for perioidc 
   
    % spread force to fluid in x direction
    ww = w(:,:,k);
    fu(i1,i2,1) = fu(i1,i2,1) + (c*F(k,1))*ww; 
    fu(i1,i2,2) = fu(i1,i2,2) + (c*F(k,2))*ww; 
   
end



end