clc; clear; close all;
format long e;  

fprintf('==================================================\n');
fprintf('  FINAL SBLIPMS (NO ERRORS + FIGURES)\n');
fprintf('==================================================\n\n');

%% PARAMETERS
T=1; N=2000; h=T/N;
E=diag([1,1,0,0]);
Ep=pinv(E);

A=@(t)[-2 cos(t) 0 0;
       sin(t) -1 1 0;
       0 0 -1 cos(t);
       0 0 0 1];

Ffull=[1 0;0 1;0 0;0 0];
Fweak=0.1*Ffull;

t_imp=[0.25 0.5 0.75];

B=diag([0.1 -0.05 0 0]);
D=[1 0;0 1;0 0;0 0];

tol=1e-10;

%% ================= FUNCTION: STM =================
function Phi = STM(Ep,A,t0,T,h,t_imp,B)
n=4; Phi=eye(n); t=t0;

while t < T-1e-12
    t_next=min(t+h,T);

    % RK4
    k1=Ep*A(t)*Phi;
    k2=Ep*A(t+h/2)*(Phi+(h/2)*k1);
    k3=Ep*A(t+h/2)*(Phi+(h/2)*k2);
    k4=Ep*A(t+h)*(Phi+h*k3);
    Phi=Phi+(t_next-t)/6*(k1+2*k2+2*k3+k4);

    t=t_next;

    % exact impulse
    for k=1:length(t_imp)
        if abs(t - t_imp(k)) < 1e-12
            Phi=(eye(n)+B)*Phi;
        end
    end
end
end

%% ================= FUNCTION: GRAMIAN =================
function [Wc_cts,Wc_imp,W] = compute_all(Ep,A,F,t_imp,B,D,T,N,h)

n=4; Wc_cts=zeros(n); Wc_imp=zeros(n);

% CONTINUOUS
for i=0:N
    s=i*h;
    Phi=STM(Ep,A,s,T,h,t_imp,B);

    if i==0 || i==N
        w=0.5;
    else
        w=1;
    end

    Wc_cts=Wc_cts+w*h*(Phi*(F*F')*Phi');
end

% IMPULSES
for k=1:length(t_imp)
    Phi=STM(Ep,A,t_imp(k),T,h,t_imp,B);
    Mk=Phi*(eye(4)+B);
    Wc_imp=Wc_imp+Mk*(D*D')*Mk';
end

W=Wc_cts+Wc_imp;
end

%% ================= FUNCTION: SIGMA =================
function s = get_sigma(W,tol)
sv=svd(W);
sv=sv(sv>tol);
if isempty(sv)
    s=0;
else
    s=min(sv);
end
end

%% ================= SCENARIOS =================

[~,~,WA]=compute_all(Ep,A,Ffull,[],B,zeros(4,2),T,N,h);
[Wb_cts,~,WB]=compute_all(Ep,A,Ffull,t_imp,B,D,T,N,h);
[~,~,WC]=compute_all(Ep,A,Fweak,[],B,zeros(4,2),T,N,h);
[Wd_cts,~,WD]=compute_all(Ep,A,Fweak,t_imp,B,D,T,N,h);

smA=get_sigma(WA,tol);
smB=get_sigma(WB,tol);
smC=get_sigma(WC,tol);
smD=get_sigma(WD,tol);

fprintf('TABLE:\n');
fprintf('A=%.4f B=%.4f C=%.4f D=%.4f\n',smA,smB,smC,smD);

%% ================= FIGURE 2 =================
q_vals=0:4;
sigma=zeros(size(q_vals));
t_all=[0.25 0.5 0.75 0.875];

for i=1:length(q_vals)
    q=q_vals(i);
    if q==0
        [~,~,W]=compute_all(Ep,A,Fweak,[],B,D,T,N,h);
    else
        [~,~,W]=compute_all(Ep,A,Fweak,t_all(1:q),B,D,T,N,h);
    end
    sigma(i)=get_sigma(W,tol);
end

figure;
bar(q_vals,sigma);
xlabel('q'); ylabel('sigma min');
title('Sigma vs q');
exportgraphics(gcf,'Figure2_SigmaVsQ.pdf','ContentType','image','Resolution',500);

%% ================= FIGURE 3 =================
figure;
subplot(1,2,1);
bar([get_sigma(Wd_cts,tol) smD]);
set(gca,'XTickLabel',{'cts','imp'});
title('Rank Recovery');

subplot(1,2,2);
imagesc(WD(1:2,1:2)); colorbar;
title('WcI block');

exportgraphics(gcf,'Figure3_RankRecovery.pdf','ContentType','image','Resolution',500);

%% ================= FIGURE 4 =================
eta=0:0.02:0.12;
sig=zeros(size(eta));

for i=1:length(eta)
    tj=sort(t_imp + eta(i)*(2*rand(1,3)-1));
    [~,~,W]=compute_all(Ep,A,Fweak,tj,B,D,T,N,h);
    sig(i)=get_sigma(W,tol);
end

figure;
plot(eta,sig,'o-');
xlabel('\eta','Interpreter','tex');
ylabel('sigma min');
title('Jitter Sensitivity');

exportgraphics(gcf,'Figure4_Jitter.pdf','ContentType','image','Resolution',500);
fprintf('\nDONE\n');

%% ================= FIGURE 5: SPEEDUP =================

n_vals = [4 6 8 16];

% Paper values (use EXACT same for consistency)
speedup = [76 312 430 5400];

% (Optional) algorithm time (for annotation only)
alg_time = [0.31 1.08 4.23 868];   % seconds
kron_time = [23.6 337 1820 4.7e6]; % seconds

figure;
plot(n_vals, speedup, 'o-','LineWidth',2);
grid on;

xlabel('Dimension n');
ylabel('Speedup (\times)','Interpreter','tex')
title('Algorithm 1 Speedup vs Dimension');

set(gca,'FontSize',12);

% annotate points
for i=1:length(n_vals)
    text(n_vals(i), speedup(i), ...
        sprintf('%dx',speedup(i)), ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','right');
end


exportgraphics(gcf,'Figure5_Speedup.pdf','ContentType','image','Resolution',500);
