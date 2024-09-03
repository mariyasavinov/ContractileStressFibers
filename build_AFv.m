function [AFv] = build_AFv(XYfiber,hq_array)

Nb = size(XYfiber,1);
tau = (XYfiber(2:end,:)-XYfiber(1:end-1,:))./vecnorm(XYfiber(2:end,:)-XYfiber(1:end-1,:),2,2);

tVk_vk = [-tau(1,1).^2./hq_array(1); -tau(2:end,1).^2./hq_array(2:end)' - tau(1:end-1,1).^2./hq_array(1:end-1)'; -tau(end,1).^2./hq_array(end)];
tVk_wk = [-tau(1,2).^2./hq_array(1); -tau(2:end,2).^2./hq_array(2:end)' - tau(1:end-1,2).^2./hq_array(1:end-1)'; -tau(end,2).^2./hq_array(end)];
tVk = zeros(2*Nb,1); tVk(1:2:end) = tVk_vk; tVk(2:2:end) = tVk_wk;

tWk_wk = [-tau(1,2).*tau(1,1)./hq_array(1)'; -tau(2:end,2).*tau(2:end,1)./hq_array(2:end)' - tau(1:end-1,2).*tau(1:end-1,1)./hq_array(1:end-1)'; - tau(end,2).*tau(end,1)./hq_array(end)'];
tWk_vk1 = tau(:,1).*tau(:,2)./hq_array';
tWk  = zeros(2*Nb,1); tWk(2:2:end) = tWk_wk; tWk(3:2:end) = tWk_vk1;

tVk1_vk1 = tau(:,1).^2./hq_array';
tVk1_wk1 = tau(:,2).^2./hq_array';
tVk1 = zeros(2*Nb,1); tVk1(3:2:end) = tVk1_vk1; tVk1(4:2:end) = tVk1_wk1;

tVkn1_wkn1 = tau(:,1).*tau(:,2)./hq_array';
tVkn1_vk = [-tau(1,1).*tau(1,2)./hq_array(1)'; - tau(2:end,1).*tau(2:end,2)./hq_array(2:end)' - tau(1:end-1,1).*tau(1:end-1,2)./hq_array(1:end-1)'; -tau(end,1).*tau(end,2)./hq_array(end)'];
tVkn1 = zeros(2*Nb,1); tVkn1(1:2:end) = tVkn1_vk; tVkn1(2:2:end) = [tVkn1_wkn1; 0];

tWkn1_vkn1 = tau(:,1).^2./hq_array';
tWkn1_wkn1 = tau(:,2).^2./hq_array';
tWkn1 = zeros(2*Nb,1); tWkn1(1:2:end)=[tWkn1_vkn1;0]; tWkn1(2:2:end) = [tWkn1_wkn1;0];

tWk1_wk1 = tau(:,2).*tau(:,1)./hq_array';
tWk1 = zeros(2*Nb,1); tWk1(4:2:end) = tWk1_wk1;

tVkn1_2_vkn1 = tau(:,1).*tau(:,2)./hq_array';
tVkn1_2 = zeros(2*Nb,1); tVkn1_2(1:2:end) = [tVkn1_2_vkn1; 0];

AFv = spdiags([tVkn1_2 tWkn1 tVkn1 tVk tWk tVk1 tWk1],[-3 -2 -1 0 1 2 3],2*Nb,2*Nb);

end
