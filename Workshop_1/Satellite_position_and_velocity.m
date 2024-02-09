function [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(time,...
    j)
%Satellite_position_and_velocity - returns ECEF Cartesian positions and
%ECEF velocities for one satellite in a constellation of 30. Simple circular
%orbits with regularly distributed satellites are modeled.
%
% Software for use by students on COMPGX04
%
% This function created 24/10/2016 by Paul Groves
% based on Satellite_positions_and_velocities.m (created 11/4/2012)
%
% Inputs:
%   time                  Current simulation time(s)
%   j                     Satellite number
% Outputs:
%   sat_r_es_e       ECEF satellite position (m) 3x1 column vector
%   sat_v_es_e       ECEF satellite velocity (m/s) 3x1 column vector
%

% Copyright 2016, Paul Groves
% License: BSD; see license.txt for details

% Constants (sone of these could be changed to inputs at a later date)
mu = 3.986004418E14; %WGS84 Earth gravitational constant (m^3 s^-2)
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
% Number of satellites in constellation
GNSS_config.no_sat = 30;
% Orbital radius of satellites (m)
GNSS_config.r_os = 2.656175E7;
% Inclination angle of satellites (deg)
GNSS_config.inclination = 55;
% Longitude offset of constellation (deg)
GNSS_config.const_delta_lambda = 0;
% Timing offset of constellation (s)
GNSS_config.const_delta_t = 0;

% Begins

% Convert inclination angle to degrees
inclination = degtorad(GNSS_config.inclination); 

% Determine orbital angular rate using (8.8)
omega_is = sqrt(mu / GNSS_config.r_os^3);

% Determine constellation time
const_time = time + GNSS_config.const_delta_t;

% (Corrected) argument of latitude
u_os_o = 2*pi*(j-1)/GNSS_config.no_sat + omega_is*const_time;
    
% Satellite position in the orbital frame from (8.14)
r_os_o = GNSS_config.r_os*[cos(u_os_o);sin(u_os_o);0];
    
% Longitude of the ascending node from (8.16)
Omega = (pi*mod(j,6)/3 + degtorad(GNSS_config.const_delta_lambda)) -...
    omega_ie*const_time;
    
% ECEF Satellite Position from (8.19)
sat_r_es_e(1,1:3) = [r_os_o(1)*cos(Omega) - r_os_o(2)*...
    cos(inclination)*sin(Omega);...
    r_os_o(1)*sin(Omega) + r_os_o(2)*cos(inclination)*cos(Omega);...
    r_os_o(2)*sin(inclination)]';
    
% Satellite velocity in the orbital frame from (8.25), noting that with
% a circular orbit r_os_o is constant and the time derivative of u_os_o
% is omega_is.
v_os_o = GNSS_config.r_os*omega_is*[-sin(u_os_o);cos(u_os_o);0];
    
% ECEF Satellite velocity from (8.26)
sat_v_es_e(1,1:3) = [v_os_o(1)*cos(Omega) - v_os_o(2)*...
    cos(inclination)*sin(Omega) + omega_ie*sat_r_es_e(2);...
    (v_os_o(1)*sin(Omega) + v_os_o(2)*cos(inclination)*cos(Omega) -...
    omega_ie*sat_r_es_e(1)); v_os_o(2)*sin(inclination)]';
   

% Ends