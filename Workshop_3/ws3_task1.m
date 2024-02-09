time_speed_heading = readmatrix("Workshop3_Speed_Heading.csv");

time = time_speed_heading(:,1);
v_ave = time_speed_heading(:,2);
heading = deg2rad(time_speed_heading(:,3));
h = 37.4; % height
time_interval = 0.5;

% Initialise varibales
v_ave_N = zeros(length(time),1);
v_ave_E = zeros(length(time),1);
latitude_DR = zeros(length(time),1);
longitude_DR = zeros(length(time),1);
R_N = zeros(length(time),1);
R_E = zeros(length(time),1);

v_ave_N_damped = zeros(length(time),1);
v_ave_E_damped = zeros(length(time),1);

% Calculating Latitude and Longitude and Damped Instantaneous DR velocites results
for i = 1: length(time)
    if i == 1 % at the first epoch (t = 0)
        v_ave_N(i) = cos(heading(i)) * v_ave(i);
        v_ave_E(i) = sin(heading(i)) * v_ave(i);
        latitude_DR(i) = deg2rad(50.424958); % mistake made: forget to transfer unit here
        longitude_DR(i) = deg2rad(-3.5957974); 
        v_ave_N_damped(i) = v_ave_N(i);
        v_ave_E_damped(i) = v_ave_E(i);
    else
        v_ave_N(i) = (1/2) * (cos(heading(i))+cos(heading(i-1))) * v_ave(i);
        v_ave_E(i) = (1/2) * (sin(heading(i))+sin(heading(i-1))) * v_ave(i);
        % use R_N from previous epoch to calculate latitude and longitude
        [R_N(i), R_E(i)] = Radii_of_curvature(latitude_DR(i-1));
        latitude_DR(i) =  latitude_DR(i-1) + (v_ave_N(i)*time_interval)/(R_N(i)+h); % result is alread in rads
        longitude_DR(i) = longitude_DR(i-1) + (v_ave_E(i)*time_interval)/((R_E(i)+h)*cos(latitude_DR(i))); % result is alread in rads
        v_ave_N_damped(i) = 1.7 * v_ave_N(i) - 0.7 * v_ave_N_damped(i-1);
        v_ave_E_damped(i) = 1.7 * v_ave_E(i) - 0.7 * v_ave_E_damped(i-1);
    end
end

latitude_DR = rad2deg(latitude_DR);
longitude_DR = rad2deg(longitude_DR);
DR_solution = [time latitude_DR longitude_DR v_ave_N_damped v_ave_E_damped];

save("DR_solution.mat", "DR_solution")







