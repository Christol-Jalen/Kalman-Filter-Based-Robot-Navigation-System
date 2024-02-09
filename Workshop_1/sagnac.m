function C = sagnac(w_ie, r_aj, c)

C = [1, w_ie*r_aj/c, 0;
    -w_ie*r_aj/c, 1, 0;
    0, 0, 1];
end