%% Task 1
%% (a)
% initialise the kalman filter state vector estimate

x_hat_plus_zero = [2447019; -5884199; -284783;...
                   184; 77; 0]; % Position + Velicity

P_plus_zero = [100 0 0 0 0 0;
               0 100 0 0 0 0;
               0 0 100 0 0 0;
               0 0 0 25 0 0;
               0 0 0 0 25 0;
               0 0 0 0 0 25;]; % error covariance matrix

%% (b)
propagation_interval = 1;
transition_matrix = [eye(3) propagation_interval*eye(3);
                     zeros(3) eye(3)];

%% (c)
S_a = 5; %acceleration power spectral density (PSD)
Q_k_minus_one = [1/3*S_a*(propagation_interval^3)*eye(3) 1/2*S_a*(propagation_interval^2)*eye(3);
             1/2*S_a*(propagation_interval^2)*eye(3) S_a*propagation_interval*eye(3)]; % system noise covariance matrix

%% (d)
x_hat_minus_one = transition_matrix*x_hat_plus_zero;

%% (e)
P_minus_one = transition_matrix * P_plus_zero * transpose(transition_matrix) + Q_k_minus_one;

%% (f)
H_one = [1 0 0 0 0 0;
       0 1 0 0 0 0;
       0 0 1 0 0 0];

%% (g)
% error standard deviation is 2.5
R_one = [6.25 0 0;
       0 6.25 0;
       0 0 6.25]; % measurement noise covariance matrix

%% (h)
K_one = P_minus_one * transpose(H_one) * (H_one * P_minus_one * transpose(H_one) + R_one)^(-1);

%% (i)
GNSS_least_squares_solutions = readmatrix("Workshop2_GNSS_Pos_ECEF.csv");
r_curve = transpose(GNSS_least_squares_solutions(1, 2:4));
z_minus_one = r_curve - x_hat_minus_one(1:3);

%% (k)
P_plus_one = (eye(6)-K_one * H_one) * P_minus_one;

%% (j)
x_hat_plus_one = x_hat_minus_one + K_one * z_minus_one;

%% (l)
[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_hat_plus_one(1:3),x_hat_plus_one(4:6));
L_b = rad2deg(L_b);
lambda_b = rad2deg(lambda_b);













