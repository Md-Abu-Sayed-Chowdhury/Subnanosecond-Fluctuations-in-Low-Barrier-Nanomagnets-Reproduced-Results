function dmdtau = llg_sun_rhs(~, m, alpha, hp, h, hs)
    m = m / norm(m);
    mx = m(1); my = m(2); mz = m(3);
  
    % Effective magnetic field: h_eff = -hp*mx*x_hat + (mz + h)*z_hat
    heff = [-hp * mx; 0; mz + h];
    
    % Spin current vector along z_hat: hs_vec = [0; 0; hs]
    hs_vec = [0; 0; hs];
    
    % Vector Torques (Sun 2000)
    t_prec = -cross(m, heff);                       % Precession torque
    t_damp = -alpha * cross(m, cross(m, heff));     % Gilbert damping torque
    t_stt  = cross(m, cross(hs_vec, m));            % Slonczewski STT (In-plane)
    t_fl   = alpha * cross(m, hs_vec);              % Field-like torque
    
    % Total derivative dm/dtau
    dmdtau = t_prec + t_damp + t_stt + t_fl;
end