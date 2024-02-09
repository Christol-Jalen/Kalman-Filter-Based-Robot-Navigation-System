load('DR_solution.mat')
GNSS_data = readmatrix("Workshop3_GNSS_Pos_Vel_NED.csv");

time_GNSS = GNSS_data(:,1);
latitude_GNSS = deg2rad(GNSS_data(:,2)); 
latitude_DR = deg2rad(DR_solution(:,2));
longitude_GNSS = deg2rad(GNSS_data(:,3));
longitude_DR = deg2rad(DR_solution(:,3));
height_GNSS = GNSS_data(:,4);
v_n_GNSS = GNSS_data(:,5);
v_n_DR = DR_solution(:,4);
v_e_GNSS = GNSS_data(:,6);
v_e_DR = DR_solution(:,5);
tau_s = 0.5; % time interval
S_DR = 0.2; % DR velocity error power spectral density
sigma_Gr = 5; % position measurement error standard deviation
sigma_Gv = 0.02; % velocity measurement error standard deviation

%Initialise variables 
x = zeros(4,1,length(time_GNSS)); % x = (v_north_error; v_east_error; latitude_error; longitude_error)*351
P= zeros(4,4,length(time_GNSS));
phi = zeros(4, 4, length(time_GNSS));
Q = zeros(4, 4, length(time_GNSS));
R = zeros(4, 4, length(time_GNSS));
K = zeros(4, 4, length(time_GNSS)); % Kalmen gain

%% constructing Kalman Filter
% for i =1
[R_N,R_E] = Radii_of_curvature(latitude_GNSS(1));

v_uncer = 0.1;
pos_uncer = 10;
P_plus_zero = [v_uncer^2 0 0 0;
            0 v_uncer^2 0 0;
            0 0 pos_uncer^2/(R_N+height_GNSS(1))^2 0;
            0 0 0 pos_uncer^2/((R_E+height_GNSS(1))^2*cos(latitude_GNSS(1))^2)];
epsilon_z_minus = zeros(4, 1, length(time_GNSS));

P(:,:,1) = P_plus_zero;
R(:,:,1) = [sigma_Gr^2/(R_N+height_GNSS(1))^2 0 0 0;
            0 sigma_Gr^2/((R_N+height_GNSS(1))^2*cos(latitude_GNSS(1))) 0 0;
            0 0 sigma_Gv^2 0;
            0 0 0 sigma_Gv^2];

H_k = [0 0 -1 0;
       0 0 0 -1;
       -1 0 0 0;
       0 -1 0 0];

K(:,:,1) = P(:,:,1)*transpose(H_k)*inv((H_k*P(:,:,1)*transpose(H_k)+R(:,:,1)));

epsilon_z_minus(:,:,1) = [rad2deg(latitude_GNSS(1)-latitude_DR(1)); 
                      rad2deg(longitude_GNSS(1)-longitude_DR(1));
                      v_n_GNSS(1)-v_n_DR(1);
                      v_e_GNSS(1)-v_e_DR(1)] - H_k*x(:,:,1);


for i = 2:length(time_GNSS)
    [R_N,R_E] = Radii_of_curvature(latitude_GNSS(i));
    % 1. transistion matrix
    phi(:,:,i-1) = [1 0 0 0;
                     0 1 0 0;
                     tau_s/(R_N+height_GNSS(i-1)) 0 1 0;
                     0 tau_s/((R_E+height_GNSS(i-1)*cos(latitude_GNSS(i-1)))) 0 1];
    % 2. system noise covariance matrix
    Q(:, :, i-1) = [S_DR*tau_s 0 (1/2)*(S_DR*tau_s^2)/(R_N+height_GNSS(i-1)) 0;
                    0 S_DR*tau_s 0 (1/2)*(S_DR*tau_s^2)/((R_E+height_GNSS(i-1))*cos(latitude_GNSS(i-1)));
                    (1/2)*(S_DR*tau_s^2)/(R_N+height_GNSS(i-1)) 0 (1/3)*(S_DR*tau_s^3)/(R_N+height_GNSS(i-1))^2 0;
                    0 (1/2)*(S_DR*tau_s^2)/((R_E+height_GNSS(i-1))*cos(latitude_GNSS(i-1))) 0 (1/3)*(S_DR*tau_s^3)/((R_E+height_GNSS(i-1))^2*cos(latitude_GNSS(i-1))^2)];
    % 3. propagate the state estimates (x_hat_minus_k)
    x(:,:,i) = phi(:,:,i-1)*x(:,:,i-1); 
    % 4. propagate the error covariance matrix (P_hat_minus_k)
    P(:,:,i) = phi(:,:,i-1)*P(:,:,i-1)*transpose(phi(:,:,i-1))+Q(:,:,i-1);
    % 5. compute measurement matrix
        % H_k is constant and is already defined outside the loop
    % 6. compute the measurement noise covariance matrix
    R(:,:,i) = [sigma_Gr^2/(R_N+height_GNSS(i))^2 0 0 0;
            0 sigma_Gr^2/((R_N+height_GNSS(i))^2*cos(latitude_GNSS(i))) 0 0;
            0 0 sigma_Gv^2 0;
            0 0 0 sigma_Gv^2];
    % 7. compute the Kalman gain matrix
    K(:,:,i) = P(:,:,i)*transpose(H_k)*inv((H_k*P(:,:,i)*transpose(H_k)+R(:,:,i)));
    % 8. formulate the measurement innovation vector
    epsilon_z_minus(:,:,i) = [rad2deg(latitude_GNSS(i)-latitude_DR(i)); 
                              rad2deg(longitude_GNSS(i)-longitude_DR(i));
                              v_n_GNSS(i)-v_n_DR(i);
                              v_e_GNSS(i)-v_e_DR(i)] - H_k*x(:,:,i);
    % 9. update the state estimates (x_hat_plus_k)
    x(:,:,i) = x(:,:,i) + K(:,:,i)*epsilon_z_minus(:,:,i);
    % 10. update the error covariance matrix (P_hat_plus_k)
    P(:,:,i) = (eye(4) - K(:,:,i)*H_k) * P(:,:,i);
end

Corrected_Solution = zeros(length(time_GNSS), 5);
Corrected_Solution(:,1) = time_GNSS;
for i = 1:length(time_GNSS)
    Corrected_Solution(i,2:5) = transpose(DR_solution(i,2:5)' - [x(3,1,i); x(4,1,i); x(1,1,i); x(2,1,i)]);
end