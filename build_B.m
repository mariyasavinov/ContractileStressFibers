function B = build_B(Nb,hq,kBfiber,FiberLoop)

Ib = speye(Nb); eb = ones(Nb,1);
C0 = spdiags([eb -2*eb eb],[-1 0 1],Nb,Nb);
if FiberLoop==1
    C0(end,1) = 1; C0(1,end) = 1;
end
C1 = circshift(C0,-1,1);
if FiberLoop~=1
    C1(end,:) = 0; %nonperiodic
end
Cn1 = circshift(C0,1,1);
if FiberLoop~=1
    Cn1(1,:) = 0; % nonperiodic
end
C = Cn1 - 2*C0 + C1;
if FiberLoop~=1
C(1,1:3) = [1, -2, 1];
C(2,1:4) = [-2, 5, -4, 1];
C(end-1,end-3:end) = [1, -4, 5, -2];
C(end,end-2:end) = [1, -2, 1];
end

B = -kBfiber*kron(speye(2),C)/hq^3;

end