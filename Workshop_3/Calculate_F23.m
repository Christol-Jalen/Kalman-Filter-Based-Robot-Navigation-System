function F23 = Calculate_F23(r_eb_e_inertial,L_b_inertial)
% Calculate_F23 - Calculates the sub-matrix F_21 in the INS/GNSS
% integration Kalman filter
%
% This function created 30/11/2016 by Paul Groves
%
% Inputs:
%   r_eb_e_inertial   Cartesian ECEF inertial position solution
%   L_b_inertial      Inertial latitude solution
%
% Outputs:
%   F23               F_21 Sub-matrix of system matrix

% Copyright 2016, Paul Groves
% License: BSD; see license.txt for details

% Constants (sone of these could be changed to inputs at a later date)

R_0 = 6378137; %WGS84 Equatorial radius in meters
e = 0.0818191908425; %WGS84 eccentricity

% Begins
geocentric_radius = R_0 / sqrt(1 - (e * sin(L_b_inertial))^2) *...
    sqrt(cos(L_b_inertial)^2 + (1 - e^2)^2 * sin(L_b_inertial)^2);

F23 = - 2 * Gravity_ECEF(r_eb_e_inertial) /...
    geocentric_radius * r_eb_e_inertial' / sqrt (r_eb_e_inertial' *...
    r_eb_e_inertial);

% Ends