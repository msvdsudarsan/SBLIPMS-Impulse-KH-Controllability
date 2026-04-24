function Mk = compute_Mk(E, Afun, Bk_list, t_imp, k, T, h)

    n=size(E,1);
    Phi=eye(n);
    t=t_imp(k);

    while t<T-1e-14
        t_next=min(t+h,T);
        idx=find(t_imp>t & t_imp<=t_next & t_imp>t_imp(k));

        if isempty(idx)
            Phi=rk4step(E,Afun,Phi,t,t_next-t);
            t=t_next;
        else
            for j=idx(:)'
                if t_imp(j)>t
                    Phi=rk4step(E,Afun,Phi,t,t_imp(j)-t);
                end
                Phi=(eye(n)+Bk_list{j})*Phi;
                t=t_imp(j);
            end
            if t_next>t
                Phi=rk4step(E,Afun,Phi,t,t_next-t);
            end
            t=t_next;
        end
    end

    Mk = Phi*(eye(n)+Bk_list{k});
end