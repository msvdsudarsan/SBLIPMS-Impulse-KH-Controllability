%% Figure1_SatelliteSV.m
clc; clear; close all;
format short g;

fprintf('==============================================\n');
fprintf('  Figure 1: Satellite SV Trajectory\n');
fprintf('  Madhyannapu & Pradheep Kumar S., 2026\n');
fprintf('==============================================\n\n');

%% PARAMETERS
T     = 1;
N_sim = 3000;
h_sim = T / N_sim;

E6   = diag([1, 1, 1, 1, 0, 0]);
Ep6  = pinv(E6);

A6fun = @(t) [ -0.8,  0.3*sin(2*pi*t),  0, 0, 0, 0;
               0.2,  -0.6,              0.1,0, 0, 0;
               0,     0.1*cos(2*pi*t), -0.5,0.2,0, 0;
               0,     0,                0.15*sin(2*pi*t), -0.4,0,0;
               0,     0,                0, 0, 0, 0;
               0,     0,                0, 0, 0, 0];

t_imp  = [0.25, 0.75];
q      = 2;

D_sat  = [1 0 0; 0 1 0; 0 0 1; 0 0 0; 0 0 0; 0 0 0];

%% SIMULATION
fprintf('Simulating X(t) over [0, T]  (N = %d steps)...\n', N_sim);

X0 = diag([0.8, 0.6, 0.4, 0.2, 0, 0]);

t_grid  = linspace(0, T, N_sim + 1);
SV_traj = zeros(6, N_sim + 1);

SV_traj(:,1) = svd(X0);

X_now = X0;
t_now = 0;

for step = 1:N_sim
    t_next = t_now + h_sim;
    in_int = find(t_imp > t_now & t_imp <= t_next);

    if isempty(in_int)
        X_now = rk4step_X(Ep6, A6fun, X_now, t_now, h_sim);
    else
        for ki = in_int(:)'
            dt_to = t_imp(ki) - t_now;
            if dt_to > 1e-14
                X_now = rk4step_X(Ep6, A6fun, X_now, t_now, dt_to);
            end
            v_k   = 0.3 * ones(3,1);
            X_now = X_now + D_sat * (v_k * v_k') * D_sat';
            t_now = t_imp(ki);
        end
        if t_next - t_now > 1e-14
            X_now = rk4step_X(Ep6, A6fun, X_now, t_now, t_next - t_now);
        end
    end

    t_now = t_next;
    SV_traj(:,step+1) = svd(X_now);
end

fprintf('  Simulation complete.\n\n');

%% FIGURE
fprintf('Generating Figure 1...\n');

T_orb_min = 90;
t_phys = t_grid * T_orb_min;

fig5 = figure('Visible','off');
set(fig5,'Units','inches','Position',[1 1 8 5]);

set(gca,'FontName','Arial','FontSize',11);
hold on;

colors = lines(6);

for ii = 1:4
    plot(t_phys, SV_traj(ii,:), '-', 'Color', colors(ii,:), 'LineWidth',2);
end

plot(t_phys, SV_traj(5,:), '--','Color',[0.6 0.6 0.6],'LineWidth',1.2);
plot(t_phys, SV_traj(6,:), '--','Color',[0.6 0.6 0.6],'LineWidth',1.2);

ymax_val = max(SV_traj(1,:))*1.15;

for ki = 1:q
    t_k_min = t_imp(ki)*T_orb_min;
    xline(t_k_min,'r--','LineWidth',2);
   text(t_k_min - 3, ymax_val*0.92, ...
    ['Thruster firing t' num2str(ki)], ...
    'HorizontalAlignment','right', ...
    'FontSize',10,'Color',[0.7 0 0]);
end

xlabel('Time (minutes, one orbital period = 90 min)','FontSize',13);
ylabel('Singular values of X(t)','FontSize',13);

title({'Figure 1: Singular Value Trajectories of Attitude Covariance X(t)',...
       'Jumps at Thruster Impulse Instants Confirm Rank Recovery'},...
       'FontSize',11);

legend({'\sigma_1','\sigma_2','\sigma_3','\sigma_4','\sigma_5, \sigma_6'}, ...
       'Interpreter','tex','Location','northeast');
grid on; box on;
xlim([0 T_orb_min]);
ylim([0 ymax_val]);

%% SAVE
drawnow;
set(fig5,'Renderer','painters');

exportgraphics(fig5,'Figure1_SatelliteSV.pdf', ...
    'ContentType','image', ...
    'Resolution',500);

fprintf('  Figure1_SatelliteSV.pdf  saved  (500 DPI)\n');
fprintf('  Copy this PDF to your LaTeX folder and compile\n\n');
fprintf('DONE\n');

%% RK4 FUNCTION
function X_out = rk4step_X(Ep, Afun, X_in, t0, dt)
    f = @(t,X) Ep*(Afun(t)*X + X*Afun(t)');
    k1 = f(t0, X_in);
    k2 = f(t0+dt/2, X_in + (dt/2)*k1);
    k3 = f(t0+dt/2, X_in + (dt/2)*k2);
    k4 = f(t0+dt,   X_in + dt*k3);
    X_out = X_in + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end