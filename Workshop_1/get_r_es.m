function r_es = get_r_es()
    Pseudo_ranges = readtable("Workshop1_Pseudo_ranges.csv");
    num_s = table2array(Pseudo_ranges(1,2:end));
    time = table2array(Pseudo_ranges(2:end,1));
    r_es = zeros(3,length(num_s), length(time));
    for i = 1: length(time)
        for j = 1:length(num_s)
            [r_es(:, j, i),v_es] = Satellite_position_and_velocity(time(i),num_s(j));
        end
    end
end
