% Duct fan airflow control - PI vs feedforward PI vs LQR with integral action
%
% Motor: FAULHABER 2342 S 024 CR, parameters from the series 2342...CR
% datasheet (edition 2026-03-23). Load: 80 mm axial impeller on the motor
% shaft; flow inferred from speed via the fan affinity law, calibrated at
% the nominal duct restriction.
%
% Assumptions: flow proportional to speed at fixed restriction (restriction
% changes enter as shaft-torque disturbances instead); nominal aero load
% linearised as viscous drag about the operating point; feedforward assumes
% a measured load torque.
%
% A. Guruprasad, July 2026

clear; clc; close all; rng(0);

%% Motor - Faulhaber 2342S024CR, 24 V
% datasheet: R 7.14 ohm, L 264 uH, kM 28.6 mNm/A, kn 334 rpm/V,
% rotor inertia 5.8 gcm^2, rated 20.5 mNm / 0.846 A / 5660 rpm
R  = 7.14;          L  = 264e-6;
Kt = 28.6e-3;       Ke = 60/(2*pi*334);     % 334 rpm/V -> V.s/rad
J_rot = 5.8e-7;
Vmax  = 24;         I_rated = 0.846;

% Kt and Ke are the same constant in SI - if they disagree, one of the
% conversions above is wrong
assert(abs(Kt - Ke)/Kt < 0.01, 'Kt/Ke mismatch - unit conversion error');
fprintf('Kt = %.4f N.m/A vs Ke = %.4f V.s/rad  (agree to %.2f%%)\n', ...
        Kt, Ke, abs(Kt - Ke)/Kt*100);

%% Load - 80 mm axial impeller, direct drive
% impeller inertia from mass and radius: 20 g at 40 mm tip radius as a
% solid disc, J = 0.5*m*r^2. estimate, and a high one - real impellers
% carry most of their mass near the hub
m_imp = 0.020;  r_imp = 0.040;
J = J_rot + 0.5*m_imp*r_imp^2;

% operating point: 4000 rpm pushing 85 m^3/h through the nominal duct
w_op  = 4000*2*pi/60;
Q_op  = 85/3600;
K_out = Q_op/w_op;              % [m^3/s per rad/s]

% damping: aero torque taken as 8 mNm at w_op (~3.4 W shaft power, sane
% for an 80 mm fan at 4000 rpm), secant-linearised, plus the 0.99 mNm
% datasheet friction torque folded in as viscous. estimates, not data.
Bd = (8.0e-3 + 0.99e-3)/w_op;

% stay inside the continuous ratings at the operating point
tau_ss = Bd*w_op;  I_ss = tau_ss/Kt;  V_ss = R*I_ss + Ke*w_op;
fprintf('operating point: %.1f mNm, %.3f A, %.2f V  (rated: 20.5 mNm, %.3f A, %d V)\n', ...
        tau_ss*1e3, I_ss, V_ss, I_rated, Vmax);
assert(tau_ss < 20.5e-3 && I_ss < I_rated && V_ss < Vmax, 'operating point exceeds ratings');

Q_target = 85;                  % [m^3/h] for display
Q_ref    = Q_op;
w_ref    = w_op;

%% Plant
% x = [i; w], u = terminal voltage, tau_d = load torque deviation from
% nominal, y = flow [m^3/s]
A = [-R/L, -Ke/L; Kt/J, -Bd/J];
B = [1/L; 0];
E = [0; -1/J];
C = [0, K_out];

fprintf('open-loop poles: %.1f, %.2f\n', eig(A));
fprintf('rank(ctrb) = %d, rank(obsv) = %d\n', rank(ctrb(A,B)), rank(obsv(A,C)));

%% LQR with integral state
% xi = int(Qref - y); plain LQR regulates to zero and can't track.
Aa = [A, zeros(2,1); -C, 0];
Ba = [B; 0];

if rank(ctrb(Aa,Ba)) < 3
    error('augmented pair not controllable - check E and C');
end

% Bryson weighting. Not optional here: i ~ O(1), w ~ O(400), xi ~ O(1e-5),
% so Q = eye(3) lets speed dominate and the integral state does nothing
% (kept below as K_naive to show exactly that). Maxima:
%   max_i   3 A     transient excursion, stall is 3.36 A
%   max_w   40      ~10% of operating speed
%   max_xi  2.8e-5  ~1 m^3/h steady error corrected in ~0.1 s
%   max_u   24      the rail
max_i = 3;  max_w = 40;  max_xi = 2.8e-5;  max_u = 24;
Qw = diag([1/max_i^2, 1/max_w^2, 1/max_xi^2]);
Rw = 1/max_u^2;

K_lqr = lqr(Aa, Ba, Qw, Rw);
fprintf('K_lqr = [%.4f %.4f %.4e]\n', K_lqr);
fprintf('LQR closed-loop poles: %.1f, %.1f, %.1f\n', eig(Aa - Ba*K_lqr));

% the naive design, kept in the comparison on purpose: equal weights leave
% the integral state effectively unweighted and its closed-loop pole lands
% at ~ -1e-4 rad/s, so it never tracks on any useful timescale
K_naive = lqr(Aa, Ba, eye(3), 1);
fprintf('naive LQR closed-loop poles: %.1f, %.2f, %.4f\n', eig(Aa - Ba*K_naive));

%% PI + feedforward
% tuned on this plant (Kp x Ki grid, scenario 1); Ki/Kp = 8 puts the PI
% zero on the mechanical pole at 8.2 rad/s. the old gains from the pump
% project don't transfer - different plant, different error scale.
Kp = 4000;
Ki = 32000;
alpha_ff = 0.86;    % knee of the sweep below: rise time is on the
                    % saturation floor there, overshoot starts ~0.90

%% Time step
% forward Euler needs dt < 2/|fastest pole|, and the binding pole is the
% closed-loop electrical one (~4e4 rad/s under LQR), not eig(A)
A_pi_cl = [A - B*Kp*C, B*Ki; -C, 0];
p_all   = [eig(A); eig(Aa - Ba*K_lqr); eig(Aa - Ba*K_naive); eig(A_pi_cl)];
dt_lim  = 2/max(abs(real(p_all)));

dt    = 2e-5;
T_sim = 4;
fprintf('dt = %.1e, Euler stability limit (closed loop) = %.2e\n', dt, dt_lim);
assert(dt < dt_lim/2, 'dt too large for fastest closed-loop pole');

t = 0:dt:T_sim;
N = numel(t);

%% Scenarios - disturbance is a shaft-torque deviation [N.m]
% 1. nominal load: setpoint step - rise, overshoot, settling
% 2. filter clogging: +4 mNm step at t = 1.5 s, about 5x the slowest
%    settling time from scenario 1, so it can't contaminate the setpoint
%    metrics - disturbance rejection
% 3. gusty inlet: continuously varying torque - persistent excitation
t_dist = [1.5 2.75];

tau_d = zeros(3,N);
tau_d(2, t > t_dist(1) & t < t_dist(2)) = 4e-3;

gust_lp = dt/0.05;      % 50 ms first-order filter on the gust noise
g = filter(gust_lp, [1, -(1-gust_lp)], randn(1,N));
tau_d(3,:) = 1.5e-3*sin(3*t) + 6e-4*(g/std(g));

names  = {'Nominal load','Filter clogging','Gusty inlet'};
ctrls  = {'PI','FF','LQR','LQR'};
labels = {'Standard PI','Feedforward PI','LQR (Bryson)','LQR (Q = I)'};
gains  = {[], [], K_lqr, K_naive};
nc     = numel(ctrls);

%% Run
Y = zeros(3,nc,N);
for s = 1:3
    for c = 1:nc
        Y(s,c,:) = simulate(tau_d(s,:), dt, N, ctrls{c}, A, B, E, C, ...
                            Kp, Ki, alpha_ff, gains{c}, Q_ref, w_ref, ...
                            R, Kt, Ke, Bd, Vmax);
    end
end

%% Metrics
% scenario 2: rise/OS/settle on the window before the disturbance fires,
% dip and recovery for the step itself. recovery band is 0.5% - the dips
% sit well inside the 2% settling band, which would read zero for everyone.
for s = 1:3
    if s == 2, i_end = find(t < t_dist(1), 1, 'last'); else, i_end = N; end

    fprintf('\n--- %s ---\n', names{s});
    fprintf('%-16s %8s %8s %10s %10s', 'ctrl','rise','OS[%]','settle','ss');
    if s == 2, fprintf(' %8s %10s','dip[%]','rec[0.5%]'); end
    fprintf('\n');

    for c = 1:nc
        y  = squeeze(Y(s,c,:))';
        yw = y(1:i_end);  tw = t(1:i_end);
        fprintf('%-16s %8.3f %8.2f %10.3f %10.2f', labels{c}, ...
                riseTime(yw,tw,Q_target), (max(yw)-Q_target)/Q_target*100, ...
                settleTime(yw,tw,Q_target,0.02), mean(y(end-round(1/dt):end)));
        if s == 2
            % dip/recovery only mean anything against a settled baseline
            pre = y(find(t < t_dist(1), 1, 'last'));
            if abs(pre - Q_target) < 0.02*Q_target
                [d,rec] = dipMetrics(y,t,Q_target,t_dist,0.005);
                fprintf(' %8.2f %10.3f', d, rec);
            else
                fprintf(' %8s %10s', '-', '-');
            end
        end
        fprintf('\n');
    end
end

%% dt convergence check - results should not move
y1 = simulate(tau_d(1,:), dt, N, 'LQR', A,B,E,C, Kp,Ki,alpha_ff,K_lqr, ...
              Q_ref,w_ref,R,Kt,Ke,Bd,Vmax);
dt2 = dt/4;  t2 = 0:dt2:T_sim;  N2 = numel(t2);
y2 = simulate(zeros(1,N2), dt2, N2, 'LQR', A,B,E,C, Kp,Ki,alpha_ff,K_lqr, ...
              Q_ref,w_ref,R,Kt,Ke,Bd,Vmax);
fprintf('\ndt check (LQR, nominal): rise %.4f @ dt=%.0e, %.4f @ dt=%.0e\n', ...
        riseTime(y1,t,Q_target), dt, riseTime(y2,t2,Q_target), dt2);

%% Plots
if ~exist('figs','dir'), mkdir('figs'); end
cols = {'b','g','r--','m-.'};

f1 = figure('Color','w','Position',[50 50 1100 850]);
if isprop(f1,'Theme'), f1.Theme = 'light'; end     % R2025b batch defaults to dark
for s = 1:3
    subplot(3,2,2*s-1); hold on; grid on;
    for c = 1:nc
        plot(t, squeeze(Y(s,c,:)), cols{c}, 'LineWidth',1.4);
    end
    yline(Q_target,'k:','Target');
    ylabel('Flow [m^3/h]'); ylim([0 110]); title(names{s});
    if s == 1, legend(labels,'Location','SouthEast'); end
    if s == 3, xlabel('Time [s]'); end

    subplot(3,2,2*s); grid on;
    area(t, tau_d(s,:)*1e3, 'FaceColor', [.85 .85 .85]);
    ylabel('\tau_d [mNm]'); ylim([-3 6]);
    if s == 3, xlabel('Time [s]'); end
end
delete(findall(f1,'Type','AxesToolbar'));
exportgraphics(f1, 'figs/comparison.png', 'Resolution', 150, 'BackgroundColor','white');

%% Zoomed views - startup and the load-step dip
f3 = figure('Color','w','Position',[50 50 1100 420]);
if isprop(f3,'Theme'), f3.Theme = 'light'; end

subplot(1,2,1); hold on; grid on;
for c = 1:nc
    plot(t, squeeze(Y(1,c,:)), cols{c}, 'LineWidth',1.4);
end
yline(Q_target,'k:','Target');
xlim([0 0.6]); ylabel('Flow [m^3/h]'); xlabel('Time [s]');
title('Startup (nominal load)');
legend(labels,'Location','SouthEast');

subplot(1,2,2); hold on; grid on;
for c = 1:nc
    plot(t, squeeze(Y(2,c,:)), cols{c}, 'LineWidth',1.4);
end
yline(Q_target,'k:','Target');
% y-limits from the three tracking controllers; the naive LQR sits near
% zero flow and just exits the frame, which is the point
zwin = t >= 1.3 & t <= 2.3;
yz = Y(2,1:3,zwin);
xlim([1.3 2.3]); ylim([min(yz(:))-0.3, max(yz(:))+0.3]);
ylabel('Flow [m^3/h]'); xlabel('Time [s]');
title('Load step at t = 1.5 s (zoom)');

delete(findall(f3,'Type','AxesToolbar'));
exportgraphics(f3, 'figs/zoom.png', 'Resolution', 150, 'BackgroundColor','white');

%% Feedforward scaling sweep
a = 0.5:0.02:1.0;
rt = zeros(size(a)); os = zeros(size(a));
for k = 1:numel(a)
    y = simulate(tau_d(1,:), dt, N, 'FF', A,B,E,C, Kp,Ki,a(k),K_lqr, ...
                 Q_ref,w_ref,R,Kt,Ke,Bd,Vmax);
    rt(k) = riseTime(y,t,Q_target);
    os(k) = (max(y)-Q_target)/Q_target*100;
end

fprintf('\nfeedforward sweep:\n%8s %8s %8s\n','alpha','rise','OS[%]');
fprintf('%8.2f %8.3f %8.2f\n', [a; rt; os]);

f2 = figure('Color','w');
if isprop(f2,'Theme'), f2.Theme = 'light'; end
yyaxis left;  plot(a, rt, 'LineWidth',1.5); ylabel('Rise time [s]');
yyaxis right; plot(a, os, 'LineWidth',1.5); ylabel('Overshoot [%]');
xlabel('\alpha'); grid on; title('Feedforward scaling tradeoff');
delete(findall(f2,'Type','AxesToolbar'));
exportgraphics(f2, 'figs/ff_sweep.png', 'Resolution', 150, 'BackgroundColor','white');


function y = simulate(tau_d, dt, N, mode, A,B,E,C, Kp,Ki,alpha_ff,K_lqr, ...
                      Q_ref, w_ref, R, Kt, Ke, Bd, Vmax)
    x = [0;0]; xi = 0; y = zeros(1,N);

    for k = 1:N-1
        tau = tau_d(k);
        e   = Q_ref - C*x;

        switch mode
            case 'PI'
                u = Kp*e + Ki*xi;
            case 'FF'
                % inverse dynamics at the reference; nominal drag Bd*w_ref
                % plus the measured disturbance torque
                V_ff = R*((Bd*w_ref + tau)/Kt) + Ke*w_ref;
                u = Kp*e + Ki*xi + alpha_ff*V_ff;
            case 'LQR'
                u = -K_lqr*[x; xi];
        end

        u_sat = min(Vmax, max(0,u));

        % conditional integration - without this the integrator charges
        % through the startup clamp and dumps it as overshoot on release
        if ~(u_sat ~= u && e*(u - u_sat) > 0)
            xi = xi + e*dt;
        end

        x = x + (A*x + B*u_sat + E*tau)*dt;
        y(k+1) = C*x;
    end

    y = y*3600;     % m^3/s -> m^3/h for display
end


function tr = riseTime(y, t, target)
    i10 = find(y >= 0.10*target, 1);
    i90 = find(y >= 0.90*target, 1);
    if isempty(i10) || isempty(i90)
        tr = NaN;
    else
        tr = t(i90) - t(i10);
    end
end


function ts = settleTime(y, t, target, tol)
    out = find(abs(y-target) > tol*target);
    if isempty(out), ts = 0; else, ts = t(out(end)); end
end


function [d, rec] = dipMetrics(y, t, target, win, tol)
    pre = y(find(t < win(1), 1, 'last'));
    idx = t >= win(1) & t <= win(2);
    d   = (pre - min(y(idx)))/pre*100;

    % onset transient only - searching past win(2) would measure the
    % release step instead
    late = find(idx & abs(y-target) > tol*target, 1, 'last');
    if isempty(late), rec = 0; else, rec = t(late) - win(1); end
end
