function [x_est,P_matrix] = Initialise_GNSS_KF
%Initialise_GNSS_KF - Initializes the GNSS EKF state estimates and error
%covariance matrix for Workshop 2
%
% This function created 30/11/2016 by Paul Groves
%
% Outputs:
%   x_est                 Kalman filter estimates:
%     Rows 1-3            estimated ECEF user position (m)
%     Rows 4-6            estimated ECEF user velocity (m/s)
%     Row 7               estimated receiver clock offset (m) 
%     Row 8               estimated receiver clock drift (m/s)
%   P_matrix              state estimation error covariance matrix

% Copyright 2016, Paul Groves
% License: BSD; see license.txt for details

% Begins

% Initialise state estimates
x_est = [  2447023.4; -5884195.9; -284784.3;...
    184.5; 76.9;  -0.1;...
    9901.1; 99.9]; 

% Initialise error covariance matrix
P_matrix =  zeros(8);
P_matrix(1,1) = 100;
P_matrix(2,2) = 100;
P_matrix(3,3) = 100;
P_matrix(4,4) = 0.01;
P_matrix(5,5) = 0.01;
P_matrix(6,6) = 0.01;
P_matrix(7,7) = 100;
P_matrix(8,8) = 0.01;

% Ends