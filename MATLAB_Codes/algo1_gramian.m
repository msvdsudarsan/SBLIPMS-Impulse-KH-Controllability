function [Wc_cts, Wc_imp, Wc_total] = ...
    algo1_gramian(E,Afun,F,Bk,Ck,Dk,t_imp,T,N)

    n=size(E,1);
    h=T/N;
    s_grid=0:h:T;

    %% ✅ CONTINUOUS (NO IMPULSES)
    Wc_cts=zeros(n);

    for i=1:length(s_grid)
        s=s_grid(i);

        Phi=eye(n);
        t=s;
        while t<T-1e-14
            dt=min(h,T-t);
            Phi=rk4step(E,Afun,Phi,t,dt);
            t=t+dt;
        end

        G = Phi*(F*F')*Phi';

        if i==1 || i==length(s_grid)
            w=0.5;
        else
            w=1;
        end

        Wc_cts = Wc_cts + w*h*G;
    end

    %% ✅ IMPULSIVE
    Wc_imp=zeros(n);
    for k=1:length(t_imp)
        Mk=compute_Mk(E,Afun,Bk,t_imp,k,T,h);
        Wc_imp=Wc_imp+Mk*(Dk{k}*Dk{k}')*Mk';
    end

    Wc_total=Wc_cts+Wc_imp;
end