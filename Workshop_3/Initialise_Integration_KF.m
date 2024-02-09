function [x_est,P_matrix] = Initialise_Integration_KF
%Initialise_Integration_KF - Initializes the Integration KF state estimates
% and error covariance matrix for Workshop 2
%
% This function created 30/11/2016 by Paul Groves
%
% Outputs:
%   x_est                 Kalman filter estimates:

%   P_matrix              state estimation error covariance matrix

% Copyright 2016, Paul Groves
% License: BSD; see license.txt for details

% Begins
deg_to_rad = 0.01745329252;
micro_g_to_meters_per_second_squared = 9.80665E-6;

% Initialise state estimates
x_est = zeros(15,1);

% Initialise error covariance matrix
P_matrix =  zeros(15);
P_matrix(1:3,1:3) = eye(3) * deg_to_rad^2;
P_matrix(4:6,4:6) = eye(3) * 0.1^2;
P_matrix(7:9,7:9) = eye(3) * 10^2;
P_matrix(10:12,10:12) = eye(3) *...
    (1000 * micro_g_to_meters_per_second_squared)^2;
P_matrix(13:15,13:15) = eye(3) * (10 * deg_to_rad / 3600)^2;

% Ends