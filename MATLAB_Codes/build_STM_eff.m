function Phi = build_STM_eff(E, Afun, Bk_list, t_imp, s, T, h)
    n = size(E,1);
    Phi = eye(n);
    t = s;

    while t < T-1e-14
        t_next = min(t+h,T);
        idx = find(t_imp>t & t_imp<=t_next & t_imp>s);

        if isempty(idx)
            Phi = rk4step(E,Afun,Phi,t,t_next-t);
            t = t_next;
        else
            for k=idx(:)'
                if t_imp(k)>t
                    Phi = rk4step(E,Afun,Phi,t,t_imp(k)-t);
                end
                Phi = (eye(n)+Bk_list{k})*Phi;
                t = t_imp(k);
            end
            if t_next>t
                Phi = rk4step(E,Afun,Phi,t,t_next-t);
            end
            t = t_next;
        end
    end
end