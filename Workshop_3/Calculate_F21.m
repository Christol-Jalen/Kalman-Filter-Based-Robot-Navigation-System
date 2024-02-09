function F21 = Calculate_F21(C_b_e_inertial,meas_f_ib_b)
% Calculate_F21 - Calculates the sub-matrix F_21 in the INS/GNSS
% integration Kalman filter
%
% This function created 30/11/2016 by Paul Groves
%
% Inputs:
%   C_b_e_inertial    Body to ECEF coordinate transformation matrix
%                     from the inertial navigation solution
%   meas_f_ib_b       Specific force measurements (in inertial navigation 
%                     solution file)
%
% Outputs:
%   F21               F_21 Sub-matrix of system matrix

% Copyright 2016, Paul Groves
% License: BSD; see license.txt for details

% Begins
F21 = -Skew_symmetric(C_b_e_inertial * meas_f_ib_b);

% Ends