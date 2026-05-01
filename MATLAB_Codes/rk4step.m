function Phi_out = rk4step(E, Afun, Phi_in, t0, dt)

    Ep = pinv(E);   % ✅ SAFE for singular E

    f  = @(t, P) Ep * (Afun(t) * P);
 
    k1 = f(t0,        Phi_in);
    k2 = f(t0+dt/2,   Phi_in + (dt/2)*k1);
    k3 = f(t0+dt/2,   Phi_in + (dt/2)*k2);
    k4 = f(t0+dt,     Phi_in +  dt   *k3);

    Phi_out = Phi_in + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end
