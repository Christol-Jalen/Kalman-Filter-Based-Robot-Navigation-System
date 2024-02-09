%latitude = -33.821075 degree = -0.5900 rad
%longitude = 151.19 degree = 2.6374 rad

%r_ea = pv_NED_to_ECEF(-0.59,2.6374,120,0);
r_ea = zeros(3,1,2); % first (3,1) current value; second (3,1) last value
satellite_j = [2, 17, 18, 22, 23, 26, 27, 28];
times = [0, 60, 120, 180, 240, 300, 360, 420, 480, 540, 600];

r_ej = zeros(3, 8);
v_ej = zeros(3, 8);
r_aj = zeros(8, 1);
r_aj_cal = zeros(8,1);
u_aj = zeros(3, 8);
C_I_e = zeros(3,3,8);
pseudo_range = readmatrix("Workshop1_Pseudo_ranges.csv");

epochs = length(pseudo_range(2:end,:));

L_b_M = zeros(length(epochs),1);
lambda_b_M = zeros(length(epochs), 1);
h_b_M = zeros(length(epochs), 1);
diff = 100;
j = 0;

for k = 1:epochs

    disp("times " + k);

    if k == 1
    
        while diff > 0.1
    
            j = j+1;
        
            disp("iteration " + j)
         
            for i = 1:8
                [r_ej(:,i),v_ej(:,i)] = Satellite_position_and_velocity(times(k),satellite_j(i));
            end
            
            
            for i = 1:8
                r_aj_cal(i) = sqrt(transpose(r_ej(:,i) - r_ea(:,:,1)) * (r_ej(:,i) - r_ea(:,:,1)));
                C_I_e(:,:,i) = sagnac(7.292115e-5, r_aj_cal(i), 299792458);
                r_aj(i) = sqrt(transpose(C_I_e(:,:,i) * r_ej(:,i)-r_ea(:,:,1)) * (C_I_e(:,:,i) * r_ej(:,i) - r_ea(:,:,1)));
                u_aj(:, i) = (C_I_e(:,:,i) * r_ej(:, i) - r_ea(:,:,1)) / r_aj(i);
            end
            
            % (e)
            offset_pre = 1;
            
            x_hat_minus = [r_ea(:,:,1); offset_pre];
            
            z_minus = zeros(8,1);
            H_e_G = ones(8,4);
            
            for i = 1:8
                z_minus(i) = pseudo_range(2,i+1) - r_aj(i) - offset_pre;
                H_e_G(i, 1) = -u_aj(1,i);
                H_e_G(i, 2) = -u_aj(2,i);
                H_e_G(i, 3) = -u_aj(3,i);
            end
            
            % (f)
            x_hat_plus = x_hat_minus + inv(transpose(H_e_G)*H_e_G)*transpose(H_e_G)*z_minus;
            % (g)
            [L_b,lambda_b,h_b,~] = pv_ECEF_to_NED(x_hat_plus(1:3,:),0);
            L_b_deg = rad2deg(L_b);
            lambda_b_deg = rad2deg(lambda_b);
            h_b;
    
            r_ea(:,:,2) = r_ea(:,:,1); % move the old result to the second (3,1)
            r_ea(:,:,1) = x_hat_plus(1:3,:); % store the new result in the first (3,1)
            r_ea
            diff = sqrt(sum(r_ea(:,:,1) - r_ea(:,:,2)).^2);
        
            [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_ea(:,:,1),0);
        
            L_b = rad2deg(L_b);
            lambda_b = rad2deg(lambda_b);
            h_b;
    
            L_b_M(k,1) = L_b;
            lambda_b_M(k,1) = lambda_b;
            h_b_M(k,1) = h_b;  
        end

    else % k != 1 ----------------------

        for i = 1:8
            [r_ej(:,i),v_ej(:,i)] = Satellite_position_and_velocity(times(k),satellite_j(i));
        end
        
        
        for i = 1:8
            r_aj_cal(i) = sqrt(transpose(r_ej(:,i) - r_ea(:,:,1)) * (r_ej(:,i) - r_ea(:,:,1)));
            C_I_e(:,:,i) = sagnac(7.292115e-5, r_aj_cal(i), 299792458);
            r_aj(i) = sqrt(transpose(C_I_e(:,:,i) * r_ej(:,i)-r_ea(:,:,1)) * (C_I_e(:,:,i) * r_ej(:,i) - r_ea(:,:,1)));
            u_aj(:, i) = (C_I_e(:,:,i) * r_ej(:, i) - r_ea(:,:,1)) / r_aj(i);
        end
        
        % (e)
        offset_pre = 1;
        
        x_hat_minus = [r_ea(:,:,1); offset_pre];
        
        z_minus = zeros(8,1);
        H_e_G = ones(8,4);
        
        for i = 1:8
            z_minus(i) = pseudo_range(2,i+1) - r_aj(i) - offset_pre;
            H_e_G(i, 1) = -u_aj(1,i);
            H_e_G(i, 2) = -u_aj(2,i);
            H_e_G(i, 3) = -u_aj(3,i);
        end
        
        % (f)
        x_hat_plus = x_hat_minus + inv(transpose(H_e_G)*H_e_G)*transpose(H_e_G)*z_minus;
        % (g)
        [L_b,lambda_b,h_b,~] = pv_ECEF_to_NED(x_hat_plus(1:3,:),0);
        L_b_deg = rad2deg(L_b);
        lambda_b_deg = rad2deg(lambda_b);
        h_b;

        r_ea(:,:,2) = r_ea(:,:,1); % move the old result to the second (3,1)
        r_ea(:,:,1) = x_hat_plus(1:3,:); % store the new result in the first (3,1)
        r_ea
        [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(r_ea(:,:,1),0);
    
        L_b = rad2deg(L_b);
        lambda_b = rad2deg(lambda_b);
        h_b;

        L_b_M(k,1) = L_b;
        lambda_b_M(k,1) = lambda_b;
        h_b_M(k,1) = h_b;  
        
    end % ------------------------------

diff = 100;
j = 0;
end
