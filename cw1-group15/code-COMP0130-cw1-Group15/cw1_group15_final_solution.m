close all;
clear;
clc;
format long g;
Define_Constants;

% read files
data_ws1_pseudo_ranges = readmatrix('Pseudo_ranges.csv');
data_ws1_pseudo_range_rates = readmatrix('Pseudo_range_rates.csv');
data_dead_reckoning = readmatrix("Dead_reckoning.csv");

satellite_list = data_ws1_pseudo_ranges(1, 2:end); % [5,6,7,9,10,11,15,30]
time_list = data_ws1_pseudo_ranges(2:end, 1); % [0:0.5:425]
[dust, sat_num] = size(satellite_list); % sat_num=8
[time_num, dust] = size(time_list); % time_num=851

%% step 1. initialise position and velocity at time=0

latitude_0 = 0; % (degree) longitude at time 0
longitude_0 = 0; % (degree) latitude at time 0
height_0 = 0; % (m) height at time 0 in NED
vel_0 = [0; 0; 0]; % (m/s) velocity at time 0 in NED

pos_offset = 0.0; % position receiver clock offset
vel_offset = 0.0; % velocity receiver clock offset

diff = 100;
epoch_time = 0;

while diff > 0.001 % threshold when terminate iteration
    % a) convert NED to ECEF
    [user_appro_pos, user_appro_vel] = pv_NED_to_ECEF(latitude_0*deg_to_rad, longitude_0*deg_to_rad, height_0, vel_0);

    % b) compute the ECEF positions and velocities of the satellites at time 0
    satellite_pos = zeros(3,sat_num);
    satellite_vel = zeros(3,sat_num);
    for i = 1:sat_num
        [satellite_pos(:,i), satellite_vel(:,i)] = Satellite_position_and_velocity(time_list(1), satellite_list(i));
    end

    % c) predict the ranges from the approximate user position to each satellite
    C = [1, 0, 0; 0, 1, 0; 0, 0, 1];
    C_I = zeros(sat_num, 3, 3);
    ranges_list = zeros(1,sat_num);
    for i = 1:sat_num
        % compute the range with the Sagnac effect compensation matrix set to the identity matrix
        ranges_list(i) = sqrt((C * satellite_pos(:,i) - user_appro_pos)' * (C * satellite_pos(:,i) - user_appro_pos));
        % compute the Sagnac effect compensation matrix
        C_I(i,:,:) = [1, omega_ie * ranges_list(i) / c, 0;
                        -omega_ie * ranges_list(i) / c, 1, 0;
                        0, 0, 1];
        % compute the ranges
        ranges_list(i) = sqrt((squeeze(C_I(i,:,:)) * satellite_pos(:,i) - user_appro_pos)' * (squeeze(C_I(i,:,:)) * satellite_pos(:,i) - user_appro_pos));
    end

    % d) compute the unit vector
    unit_vector = zeros(3,sat_num);
    for i = 1:sat_num
        unit_vector(:,i) = (squeeze(C_I(i,:,:)) * satellite_pos(:,i) - user_appro_pos) / ranges_list(i);
    end

    % predict the range rates from the approximate user velocity to each satellite
    range_rates_list = zeros(1,sat_num);
    for i = 1:sat_num
        range_rates_list(i) = unit_vector(:,i)' * ((squeeze(C_I(i,:,:)) * (satellite_vel(:,i) + Omega_ie * satellite_pos(:,i)) - (user_appro_vel + Omega_ie * user_appro_pos)));
    end

    % e) formulate the predicted state vector, measurement innovation vector and measurement matrix
    x_hat_minus = zeros(4,1); % predicted state vector for position
    x_hat_minus(1:3, 1) = user_appro_pos;
    x_hat_minus(4,1) = pos_offset;

    delta_z_minus = zeros(sat_num,1); % measurement innovation vector for position
    rho_a_j = data_ws1_pseudo_ranges(2, 2:end); % measured range from satellite to the user antenna
    for i = 1:sat_num
        delta_z_minus(i) = rho_a_j(i) - ranges_list(i) - pos_offset;
    end

    H = zeros(sat_num,4); % measurement matrix
    for i = 1:sat_num
        H(i, 1:3) = -unit_vector(:,i);
        H(i, 4) = 1;
    end

    x_dot_hat_minus = zeros(4,1); % predicted state vector for velocity
    x_dot_hat_minus(1:3, 1) = user_appro_vel;
    x_dot_hat_minus(4,1) = vel_offset;

    delta_z_dot_minus = zeros(sat_num,1); % measurement innovation vector for velocity
    rho_dot_a_j = data_ws1_pseudo_range_rates(2, 2:end); % measured range rates from satellite to the user antenna
    for i = 1:sat_num
        delta_z_dot_minus(i) = rho_dot_a_j(i) - range_rates_list(i) - vel_offset;
    end

    % f) calculate the results
    pos_cal = x_hat_minus + (H' * H) \ H' * delta_z_minus;
    vel_cal = x_dot_hat_minus + (H' * H) \H' * delta_z_dot_minus;

    % g) change to NED frame and update clock offset
    [lat_exact, lon_exact, height_exact, vel_exact] = pv_ECEF_to_NED(pos_cal(1:3), vel_cal(1:3));
    lat_exact = lat_exact * rad_to_deg;
    lon_exact = lon_exact * rad_to_deg;
    pos_offset = pos_cal(4);
    vel_offset = vel_cal(4);

    % h) prepare for next iteration
    diff = sqrt((lat_exact-latitude_0)^2 + (lon_exact-longitude_0)^2 + (height_exact-height_0)^2 + sum((vel_exact-vel_0).^2));
    latitude_0 = lat_exact;
    longitude_0 = lon_exact;
    height_0 = height_exact;
    vel_0 = vel_exact;
    epoch_time = epoch_time + 1;
end

disp('t=0:');
disp(latitude_0);
disp(longitude_0);
disp(height_0);
disp(vel_0);
disp(pos_offset);
disp(vel_offset);
vel_0_NED = vel_0;

% final initial position and velocity for step 2 and step 4
[pos_0, vel_0] = pv_NED_to_ECEF(latitude_0*deg_to_rad, longitude_0*deg_to_rad, height_0, vel_0);

%% step 2. kalman filter based GNSS solution with outlier detection

% constants
tau_s = 0.5; % (s), propagation interval
S_ae = 0.01; % (m^2s^-3), acceleration PSD
S_cphia = 0.01; % (m^2s^-3), clock phase PSD
S_cfa = 0.04; % (m^2s^-3), clock frequency PSD
sigma_pho = 10; % (m), range measurements standard deviation error
sigma_r = 0.05; % (m/s), range rate measurements standard deviation error

sigma_rho = 5; % (m) measurement standard deviation error for outlier detection
T = 6; % outlier detection threshold

GNSS_result = zeros(time_num, 7); % [time, latitude, longitude, height, vel_n, vel_e, vel_d] * 851
GNSS_result(1, 1) = 0;
GNSS_result(1, 2) = latitude_0;
GNSS_result(1, 3) = longitude_0;
GNSS_result(1, 4) = height_0;
GNSS_result(1, 5:7) = vel_0';

% a) initialise the state vector and covariance matrix
x_pre_plus = [pos_0; vel_0; pos_offset; vel_offset];
P_pre_plus = [100 * eye(3), zeros(3,5);
              zeros(3,3), 0.01 * eye(3), zeros(3,2);
              zeros(1,6), 100, zeros(1,1);
              zeros(1,7), 0.01];

% b) compute the transition matrix
trans_mat = [eye(3), tau_s * eye(3), zeros(3,1), zeros(3,1);
             zeros(3,3), eye(3), zeros(3,1), zeros(3,1);
             zeros(1,3), zeros(1,3), 1, tau_s;
             zeros(1,3), zeros(1,3), 0, 1];

% c) compute the system noise covariance matrix
noise_mat = [1/3 * S_ae * tau_s^3 * eye(3), 1/2 * S_ae * tau_s^2 * eye(3), zeros(3,1), zeros(3,1);
             1/2 * S_ae * tau_s^2 * eye(3), S_ae * tau_s * eye(3), zeros(3,1), zeros(3,1);
             zeros(1,3), zeros(1,3), S_cphia * tau_s + 1/3 * S_cfa * tau_s^3, 1/2 * S_cfa * tau_s^2;
             zeros(1,3), zeros(1,3), 1/2 * S_cfa * tau_s^2, S_cfa * tau_s];

for iter = 2:time_num
    % d) use the transition matrix to propagate the state estimates
    x_cur_minus = trans_mat * x_pre_plus;

    % e) use this to propagate the error covariance matrix
    P_cur_minus = trans_mat * P_pre_plus * trans_mat' + noise_mat;

    % outlier detection
    is_outlier = 1;
    outlier_num = 0;
    residual_pre = 0;

    pos_j = zeros(3, sat_num); % Cartesian ECEF position of satellites at current time
    vel_j = zeros(3, sat_num); % Cartesian ECEF velocity of satellites at current time
    for j = 1:sat_num
        [pos_j(:,j), vel_j(:,j)] = Satellite_position_and_velocity(time_list(iter), satellite_list(j));
    end
    
    rho_a_j = data_ws1_pseudo_ranges(iter+1, 2:end); % measured range from satellites to the user antenna
    rho_dot_a_j = data_ws1_pseudo_range_rates(iter+1, 2:end); % measured range rates from satellites to the user antenna

    while is_outlier == 1 % loop until no outlier satellite is detected
        is_outlier = 0;

        % f) predict the ranges from the approximate user position to each satellite
        C = [1, 0, 0; 0, 1, 0; 0, 0, 1];
        C_eI = zeros(sat_num-outlier_num, 3, 3);

        r_aj_minus = zeros(1,sat_num-outlier_num);
        for j = 1:sat_num-outlier_num
            r_aj_minus(j) = sqrt((C * pos_j(:,j) - x_cur_minus(1:3))' * (C * pos_j(:,j) - x_cur_minus(1:3)));
            C_eI(j, :, :) = [1, omega_ie * r_aj_minus(j) / c, 0;
                            -omega_ie * r_aj_minus(j) / c, 1, 0;
                            0, 0, 1];
            r_aj_minus(j) = sqrt((squeeze(C_eI(j,:,:)) * pos_j(:,j) - x_cur_minus(1:3))' * (squeeze(C_eI(j,:,:)) * pos_j(:,j) - x_cur_minus(1:3)));
        end

        % g) compute the line-of-sight unit vector
        u_aj_e = zeros(3, sat_num-outlier_num);
        for j = 1:sat_num-outlier_num
            u_aj_e(:, j) = (squeeze(C_eI(j,:,:)) * pos_j(:,j) - x_cur_minus(1:3)) / r_aj_minus(j);
        end

        % h) predict the range rates
        r_aj_minus_dot = zeros(sat_num-outlier_num);
        for j = 1:sat_num-outlier_num
            r_aj_minus_dot(j) = u_aj_e(:,j)' * (squeeze(C_eI(j,:,:)) * (vel_j(:,j) + Omega_ie * pos_j(:,j)) - (x_cur_minus(4:6) + Omega_ie * x_cur_minus(1:3)));
        end

        % i) compute the measurement matrix
        H_cur = zeros((sat_num-outlier_num) * 2, 8);
        for j = 1:sat_num-outlier_num
            H_cur(j,1:3) = -u_aj_e(:, j)';
            H_cur(j+sat_num-outlier_num, 4:6) = -u_aj_e(:,j)';
            H_cur(j,7) = 1;
            H_cur(j+sat_num-outlier_num,8) = 1;
        end

        % j) compute the measurement noise covariance matrix
        R_cur = [sigma_pho^2 * eye(sat_num-outlier_num), zeros(sat_num-outlier_num, sat_num-outlier_num);
                zeros(sat_num-outlier_num, sat_num-outlier_num), sigma_r^2 * eye(sat_num-outlier_num)];

        % k) compute the kalman gain matrix
        kalman_gain = P_cur_minus * H_cur' * inv(H_cur * P_cur_minus * H_cur' + R_cur);

        % l) formulate the measurement innovation vector
        delta_z_minus = zeros((sat_num-outlier_num) * 2, 1);
        for j=1:sat_num-outlier_num
            delta_z_minus(j) = rho_a_j(j) - r_aj_minus(j) - x_cur_minus(7);
            delta_z_minus(j+sat_num-outlier_num) = rho_dot_a_j(j) - r_aj_minus_dot(j) - x_cur_minus(8);
        end
        
        % outlier detection part
        % compute the residuals vector
        v = (H_cur * inv(H_cur' * H_cur) * H_cur' - eye((sat_num-outlier_num) * 2)) * delta_z_minus;
        % compute the residuals covariance matrix
        C_v = (eye((sat_num-outlier_num) * 2) - H_cur * inv(H_cur' * H_cur) * H_cur') * sigma_rho^2;
        for k = 1:length(v)
            residual = abs(v(k)) - sqrt(C_v(k,k))*T;
            % locate the time and satellite with largest residual
            if residual > residual_pre
                idx = k;
                residual_pre = residual;
                is_outlier = 1;
                disp(['outlier occurred to satellite ' num2str(satellite_list(idx)) ' at time ' num2str(time_list(iter))]);
                % pause(0.25);
            end
        end
        if is_outlier == 1 % remove the data of the outlier satellite
            outlier_num = outlier_num + 1;
            pos_j(:,idx) = [];
            vel_j(:,idx) = [];
            rho_a_j(idx) = [];
            rho_dot_a_j(idx) = [];
        end

    end % while is_outlier == 1

    % m) update the state estimate
    x_cur_plus = x_cur_minus + kalman_gain * delta_z_minus;

    % n) update the error covariance matrix
    P_cur_plus = (eye(8) - kalman_gain * H_cur) * P_cur_minus;

    % o) convert to NED
    [lat_rad, long_rad, height, vel] = pv_ECEF_to_NED(x_cur_plus(1:3), x_cur_plus(4:6));
    lat = rad_to_deg * lat_rad;
    long = rad_to_deg * long_rad;

    % prepare for next loop
    x_pre_plus = x_cur_plus;
    P_pre_plus = P_cur_plus;
    % store the result
    GNSS_result(iter, 1) = time_list(iter);
    GNSS_result(iter, 2) = lat;
    GNSS_result(iter, 3) = long;
    GNSS_result(iter, 4) = height;
    GNSS_result(iter, 5:7) = vel';
end

% Plot the result
figure;
quiver(GNSS_result(:,3), GNSS_result(:,2), GNSS_result(:,6), GNSS_result(:,5), 0.5);
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('GNSS Velocity Vector Field');
xlim([min(GNSS_result(:,3)) max(GNSS_result(:,3))]);
ylim([min(GNSS_result(:,2)) max(GNSS_result(:,2))]);
grid on;

%% step 3. kalman filter based heading correction

% constants
sigma_gb = 1; % Gyroscope bias std deviation
sigma_gn = rad2deg(1e-4); % Gyroscope ramdom noise std deviation
sigma_mag = 4; % Magnetic compass noise error std deviation 
S_rg = (3e-6)*rad_to_deg^2; % Gyroscope noise PSD
S_bgd = (2e-7)*rad_to_deg^2; % Gyroscope bias PSD

% extract data from tables
Omega = data_dead_reckoning(:,6); % Angular rate in radius gyroscpoe
Omega = rad2deg(Omega); % Angular rate in degree
Psi_mag = data_dead_reckoning(:,7); % Heading measurement from magnetic compass

% initialization
x = zeros(2,1); % Initial state, Gyro-derived heading error and gyro bias
Psi_gyro = zeros(length(time_list),1); % Gyro-derived heading
Psi_gyro(1) = Psi_mag(1); % Use the initial heading from the compass
Psi_gyro_c = Psi_gyro; % Corrected gyroscope heading

% state estimation error covariance matrix
P = [sigma_gn^2,0;
     0,sigma_gb^2];

% state transition matrix
Phi = [1, tau_s;
       0, 1];

% System noise covariance
Q = [S_rg*tau_s+(S_bgd*tau_s^3)/3, 0.5*S_bgd*tau_s^2;
    0.5*S_bgd*tau_s^2, S_bgd*tau_s]; 

% Measurement matrix
H = [-1, 0]; 

% Measurement noise covarince
R = sigma_mag^2; 

% kalman filter iteration
for i = 2:length(time_list)
    Psi_gyro(i) = Psi_gyro(i-1) + Omega(i) * tau_s;

    % Propagate state estimate
    x_p = Phi*x;

    % Propagate error covariance
    P_p = Phi*P*Phi' + Q;

    % Kalman gain matrix
    K = P_p*H'/(H*P_p*H'+R);

    % Measurement innovation vector
    epsilon_z_pre = Psi_mag(i) - Psi_gyro(i) - H*x_p;

    % Update the state estimates
    x_u = x_p + K*epsilon_z_pre;

    % Update the error covariance matrix
    P_u = (1 - K*H)*P_p;

    % Correct the gyroscope solution with the Kalman filter
    Psi_gyro_c(i,:) = Psi_gyro(i,:) - x_u(1);

    % Update states and covariance for iteration
    x = x_u;
    P = P_u;
end

Psi_gyro_c_deg = Psi_gyro_c; % Store heading in degree for plotting

%% step 4. DR solution

% constants
d_lf = 0.4; % Distance between left and right wheel
d_am = 0.05; % Distance between center of mass and the antenna

L_k0 = deg2rad(latitude_0);
Lamda_k0 = deg2rad(longitude_0);
h = GNSS_result(:,4);

Psi_gyro_c = deg2rad(Psi_gyro_c); % Gyroscope heading corrected by magnetic compass

v_f_a = (data_dead_reckoning(:, 4) + data_dead_reckoning(:, 5))/2; % Average speed of driving wheels
yaw_rate = (data_dead_reckoning(:, 4) - data_dead_reckoning(:, 5))/d_lf; % Yaw rate
v_ant_yaw = yaw_rate * d_am; % Lateral velocity at the antenna due to yaw
v_ant_a = sqrt(v_f_a.^2 + v_ant_yaw.^2); % Average speed at the antenna

% initialization
v_a = zeros(2,1); % NED speed
L_k = zeros(length(time_list),1); % Latitude
L_k(1) = L_k0;
Lamda_k = zeros(length(time_list),1); % Longitude
Lamda_k(1) = Lamda_k0;
v_i = zeros(length(time_list),2); % Instantaneous velocity
v_i(1,1) = vel_0(1);
v_i(1,2) = vel_0(2);

% dr iteration
for i = 2:length(time_list)
    V_a2NED = [cos(Psi_gyro_c(i)) + cos(Psi_gyro_c(i-1));
                sin(Psi_gyro_c(i)) + sin(Psi_gyro_c(i-1))];
    % compute average velocity
    v_a = 0.5 * V_a2NED * v_ant_a(i);
    [R_N,R_E] = Radii_of_curvature(L_k(i-1));
    %  compute latitude and longitude
    L_k(i) = L_k(i-1) + v_a(1)*tau_s/(R_N + h(i-1));
    Lamda_k(i) = Lamda_k(i-1) + v_a(2)*tau_s/((R_E + h(i-1))*cos(L_k(i)));
    % update instantaneous speed
    v_i(i,1) = 1.7*v_a(1) - 0.7*v_i(i-1,1);
    v_i(i,2) = 1.7*v_a(2) - 0.7*v_i(i-1,2);
end

% Convert latitude and longitude to degrees
L_k = rad2deg(L_k);
Lamda_k = rad2deg(Lamda_k);

DR_result = [time_list, L_k, Lamda_k, v_i, rad2deg(Psi_gyro_c)];

% Plot the result
figure;
quiver(Lamda_k, L_k, v_i(:,2), v_i(:,1), 0.5);
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('DR Velocity Vector Field');
xlim([min(Lamda_k) max(Lamda_k)]);
ylim([min(L_k) max(L_k)]);
grid on;

%% step 5. Integrated DR/GNSS solution based on Kalman Filter
 
% constants
S_DR = 0.01; % (m^2s^-3) DR velocity error power spectral density
sigma_Gr = 10;
sigma_Gv = 0.05;
sigma_v = 0.1; % (m/s) initial velocity uncertainty
sigma_r = 10; % (m) initial position uncertainty

L_pre = GNSS_result(1,2) * deg_to_rad;
[R_N, R_E] = Radii_of_curvature(L_pre);
h_pre = GNSS_result(1,4);

% store integration result
Final_result = zeros(time_num, 6); % [time, latitude, longitude, speed_n, speed_e, heading]
Final_result(1,1) = 0;
Final_result(1,2) = GNSS_result(1,2);
Final_result(1,3) = GNSS_result(1,3);
Final_result(1,4) = GNSS_result(1,5);
Final_result(1,5) = GNSS_result(1,6);
Final_result(1,6) = DR_result(1,6);

% state vector
x_pre_plus = [0; 0; 0; 0];

% state estimation error covariance matrix
P_pre_plus = [sigma_v^2, 0, 0, 0;
              0, sigma_v^2, 0, 0;
              0, 0, sigma_r^2 / (R_N + h_pre)^2, 0;
              0, 0, 0, sigma_r^2 / ((R_E + h_pre)^2 * cos(L_pre)^2)];

for i = 2:time_num
    % 1. compute the transition matrix
    trans_mat = [1, 0, 0, 0;
                 0, 1, 0, 0;
                 tau_s / (R_N + h_pre), 0, 1, 0;
                 0, tau_s / ((R_E + h_pre) * cos(L_pre)), 0, 1];

    % 2. compute the system noise covariance matrix
    Q_pre = [S_DR * tau_s, 0, 1/2 * S_DR * tau_s^2 / (R_N + h_pre), 0;
             0, S_DR * tau_s, 0, 1/2 * S_DR * tau_s^2 / ((R_E + h_pre) * cos(L_pre));
              1/2 * S_DR * tau_s^2 / (R_N + h_pre), 0, 1/3 * S_DR * tau_s^3 / (R_N + h_pre)^2, 0;
              0, 1/2 * S_DR * tau_s^2 / ((R_E + h_pre) * cos(L_pre)), 0, 1/3 * S_DR * tau_s^3 / ((R_E + h_pre)^2 * cos(L_pre)^2)];

    % 3. propagate the state estimates
    x_cur_minus = trans_mat * x_pre_plus;

    % 4. propagate the error covariance matrix
    P_cur_minus = trans_mat * P_pre_plus * trans_mat' + Q_pre;

    % 5. compute the measurement matrix
    H_cur = [0, 0, -1, 0;
             0, 0, 0, -1;
             -1, 0, 0, 0;
             0, -1, 0, 0];

    % 6. compute the measurement noise covariance matrix
    L_cur = GNSS_result(i,2) * deg_to_rad;
    [R_N, R_E] = Radii_of_curvature(L_cur);
    h_cur = GNSS_result(i,4);

    R_cur = [sigma_Gr^2 / (R_N + h_cur)^2, 0, 0, 0;
             0, sigma_Gr^2 / ((R_E + h_cur)^2 * cos(L_cur)^2), 0, 0;
             0, 0, sigma_Gv^2, 0;
             0, 0, 0, sigma_Gv^2];

    % 7. compute the Kalman gain matrix
    K_cur = P_cur_minus * H_cur' * inv(H_cur * P_cur_minus * H_cur' + R_cur);

    % 8. formulate the measurement innovation vector
    delta_z_cur_minus = [(GNSS_result(i,2) - DR_result(i,2)) * deg_to_rad;
                         (GNSS_result(i,3) - DR_result(i,3)) * deg_to_rad;
                         GNSS_result(i,5) - DR_result(i,4);
                         GNSS_result(i,6) - DR_result(i,5)] - H_cur * x_cur_minus;

    % 9. update the state estimates
    x_cur_plus = x_cur_minus + K_cur * delta_z_cur_minus;

    % 10. update the error covariance matrix
    P_cur_plus = (eye(4) - K_cur * H_cur) * P_cur_minus;

    % Final: use the Kalman filter estimates to correct the DR solution
    Final_result(i,1) = GNSS_result(i,1);
    Final_result(i,2) = (DR_result(i,2) * deg_to_rad - x_cur_plus(3)) * rad_to_deg;
    Final_result(i,3) = (DR_result(i,3) * deg_to_rad - x_cur_plus(4)) * rad_to_deg;
    Final_result(i,4) = DR_result(i,4) - x_cur_plus(1);
    Final_result(i,5) = DR_result(i,5) - x_cur_plus(2);
    Final_result(i,6) = DR_result(i,6);
    
    % prepare for next loop
    x_pre_plus = x_cur_plus;
    P_pre_plus = P_cur_plus;
    L_pre = L_cur;
    h_pre = h_cur;
end

% Plot the GNSS/DR result
figure;
quiver(Final_result(:,3), Final_result(:,2), Final_result(:,5), Final_result(:,4), 0.5);
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('GNSS/DR Velocity Vector Field');
xlim([min(Final_result(:,3)) max(Final_result(:,3))]);
ylim([min(Final_result(:,2)) max(Final_result(:,2))]);
grid on;

% Plot the 3 results
figure('Position', [100, 100, 800, 600]);
hold on;
quiver(GNSS_result(:,3), GNSS_result(:,2), GNSS_result(:,6), GNSS_result(:,5), 'Color', 'r', 'LineWidth', 0.25);
quiver(DR_result(:,3), DR_result(:,2), DR_result(:,5), DR_result(:,4), 'Color', 'g', 'LineWidth', 0.25);
quiver(Final_result(:,3), Final_result(:,2), Final_result(:,5), Final_result(:,4), 'Color', 'b', 'LineWidth', 1.5);
legend('GNSS', 'DR', 'GNSS/DR');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('GNSS/DR Velocity Vector Field');
grid on;

% Plot the heading
figure;
plot(time_list,Psi_gyro, 'g--','LineWidth', 1);
hold on;
plot(time_list,Psi_mag, 'r--','LineWidth', 2);
hold on;
plot(time_list,Psi_gyro_c_deg, 'b-','LineWidth', 1);
hold on;
xlabel('time (s)');
ylabel('Headings (degree)');
legend('Gyro-derived Heading','Magnetic Compass Measured Heading', 'Corrected Gyro Heading');
title('Gyro-Magnetometer Integration');
grid on;

% Plot the speed of the wheels
figure;
plot(time_list,data_dead_reckoning(:, 2), 'g-','LineWidth', 1);
hold on;
plot(time_list,data_dead_reckoning(:, 3), 'r-','LineWidth', 1);
hold on;
plot(time_list,data_dead_reckoning(:, 4), 'b-','LineWidth', 1);
hold on;
plot(time_list,data_dead_reckoning(:, 5), 'y-','LineWidth', 1);
hold on;
xlabel('time (s)');
ylabel('Speed of wheels (m/s)');
legend('Wheel front left','Wheel front right', 'Wheel rear left', 'Wheel rear right');
title('Speed of the wheels');
grid on;
hold off;

% Save the result
fileID = fopen('Result_Group15.csv', 'w');
fprintf(fileID, '%.1f,%.8f,%.8f,%.8f,%.8f,%.8f\n', Final_result.');
fclose(fileID);