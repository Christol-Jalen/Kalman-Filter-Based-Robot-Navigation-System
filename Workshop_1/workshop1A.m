%latitude = -33.821075 degree = -0.5900 rad
%longitude = 151.19 degree = 2.6374 rad

[r_ea, v_ea] = pv_NED_to_ECEF(deg2rad(-33.821075),deg2rad(151.188496),120,0);
satellite_j = [2, 17, 18, 22, 23, 26, 27, 28];

r_ej = get_r_es(); %8*3*11
r_ej_zero = r_ej(:, :, 1); %8*3
v_ej = zeros(3, 8);
r_aj = zeros(8, 1);
u_aj = zeros(3, 8);
C_i_e = zeros(3,3,8);


for i = 1:8 % For each satellite
    [r_ej_zero(:,i),v_ej(:,i)] = Satellite_position_and_velocity(0,satellite_j(i));
end

%r_ej = transpose (r_ej);
%v_ej = transpose (v_ej);


for i = 1:8 % For each satellite
    r_aj(i,:) = sqrt(transpose(r_ej_zero(:,i) - r_ea) * (r_ej_zero(:,i) - r_ea));
    C_i_e(:,:,i) = sagnac(7.292115e-5, r_aj(i), 299792458);
    r_aj(i,:) = sqrt(transpose(C_i_e(:,:,i) * r_ej_zero(:,i)-r_ea) * (C_i_e(:,:,i) * r_ej_zero(:,i) - r_ea));
    u_aj(:, i) = (C_i_e(:,:,i) * r_ej_zero(:, i) - r_ea) / r_aj(i);
end

% (e)
pseudo_range = readmatrix("Workshop1_Pseudo_ranges.csv");
offset_pre = 0;

x_hat_minus = [r_ea; offset_pre];
%x_hat_minus = r_ea;

epsilon_z_minus = zeros(8,1);
H_e_G = ones(8,4);

for i = 1:8
    epsilon_z_minus(i) = pseudo_range(2,i+1) - r_aj(i) - offset_pre;
    H_e_G(i, 1) = -u_aj(1,i);
    H_e_G(i, 2) = -u_aj(2,i);
    H_e_G(i, 3) = -u_aj(3,i);
end

% (f)
x_hat_plus = x_hat_minus + inv(transpose(H_e_G)*H_e_G)*transpose(H_e_G)*epsilon_z_minus;

% (g)
[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_hat_plus(1:3,:),0);
L_b = rad2deg(L_b);
lambda_b = rad2deg(lambda_b);



