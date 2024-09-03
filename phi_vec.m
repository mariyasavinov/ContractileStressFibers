function phi_r = phi_vec(r,indx)

phi_r = zeros(4,length(r),4);
q = sqrt(1+4*r.*(1-r));
if indx == 1
phi_r(4,:,:) = repmat((1+2*r-q)/8,[1,1,4]);
phi_r(3,:,:) = repmat((1+2*r+q)/8,[1,1,4]);
phi_r(2,:,:) = repmat((3-2*r+q)/8,[1,1,4]);
phi_r(1,:,:) = repmat((3-2*r-q)/8,[1,1,4]);
else
phi_r(:,:,4) = repmat((1+2*r-q)'/8,[4,1]);
phi_r(:,:,3) = repmat((1+2*r+q)'/8,[4,1]);
phi_r(:,:,2) = repmat((3-2*r+q)'/8,[4,1]);
phi_r(:,:,1) = repmat((3-2*r-q)'/8,[4,1]);
end

end