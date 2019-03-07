
tot_vol = 1.9e-3_rprec/60._rprec

!! for uniform distribution v_i = v/15

do ip = 1,npcon
    q_src(ipcon) = v_i*6/pi/diameter(ip)**3
end

