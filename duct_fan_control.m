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
% A. Guruprasad, June 2026

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
% NB: fprintf with %f drops imaginary parts silently - print poles through
% poleStr so a complex pair can't masquerade as a repeated real pole
fprintf('LQR closed-loop poles: %s\n', poleStr(eig(Aa - Ba*K_lqr)));

% the naive design, kept in the comparison on purpose: equal weights leave
% the integral state effectively unweighted and its closed-loop pole lands
% at -5.6e-5 rad/s, so it never tracks on any useful timescale
K_naive = lqr(Aa, Ba, eye(3), 1);
fprintf('naive LQR closed-loop poles: %s\n', poleStr(eig(Aa - Ba*K_naive)));

%% Luenberger observer - estimates current from the flow measurement
% (A,C) is observable (checked above), so full state feedback does not
% require a current sensor. Observer poles chosen faster than the
% controller's dominant pair but not so split that the slow mode lags
% through startup saturation - see the pole sweep below.
L_obs = place(A', C', [-2000, -20000])';
fprintf('L_obs = [%.4e %.4e]\n', L_obs);
fprintf('observer poles: %s\n', poleStr(eig(A - L_obs*C)));

%% Pole placement with static prefilter - the no-integrator strawman
% Bass-Gura placement at the LQR's pole real parts, u = -K_pp*x + G_pp*ref
% with G_pp a static gain inverting the nominal DC gain. Tracks exactly on
% the nominal plant; has no mechanism to notice the plant changed, so any
% load deviation becomes a permanent offset - the case for the integral
% state, shown rather than asserted.
K_pp = place(A, B, [-40612, -78.89]);
G_pp = -1/(C/(A - B*K_pp)*B);
fprintf('K_pp = [%.4f %.4f], G_pp = %.4e\n', K_pp, G_pp);

%% PI design - analytic, from the reduced first-order plant
% Timescale separation (electrical pole 3300x faster than mechanical)
% justifies di/dt -> 0, giving P(s) = kdc/(tau*s + 1) with the back-EMF
% damping Kt*Ke/R folded into the effective friction. PI then yields an
% exactly second-order closed loop:
%     Kp = (2*zeta*wn*tau - 1)/kdc
%     Ki = wn^2*tau/kdc
% Note this ignores input saturation, which does bind during startup;
% the design holds because saturation clears early.
B_eff = Bd + Kt*Ke/R;          % 1.360e-4 - 84% of it is electrical
tau_m = J/B_eff;               % 0.1219 s  (1/tau_m = 8.20 rad/s)
kdc   = K_out*Kt/(R*B_eff);    % 1.660e-3 m^3/s per V

zeta_d = 1.50;  wn_d = 20.9;   % deliberately overdamped baseline
Kp = (2*zeta_d*wn_d*tau_m - 1)/kdc;
Ki = wn_d^2*tau_m/kdc;
fprintf('PI design: zeta = %.2f, wn = %.1f -> Kp = %.0f, Ki = %.0f\n', ...
        zeta_d, wn_d, Kp, Ki);

alpha_ff = 0.86;    % knee of the sweep below: rise time is on the
                    % saturation floor there, overshoot starts ~0.90

%% Time step
% forward Euler needs dt < 2/|fastest pole|, and the binding pole is the
% closed-loop electrical one (~4e4 rad/s under LQR), not eig(A)
A_pi_cl = [A - B*Kp*C, B*Ki; -C, 0];
p_all   = [eig(A); eig(Aa - Ba*K_lqr); eig(Aa - Ba*K_naive); eig(A_pi_cl); ...
           eig(A - L_obs*C); eig(A - B*K_pp)];
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
ctrls  = {'PI','FF','LQR','LQR','OBS','PP'};
labels = {'Standard PI','Feedforward PI','LQR (Bryson)','LQR (Q = I)', ...
          'LQR + observer','Pole place + prefilter'};
gains  = {[], [], K_lqr, K_naive, K_lqr, K_pp};
nc     = numel(ctrls);

%% Run
Y = zeros(3,nc,N);
for s = 1:3
    for c = 1:nc
        Y(s,c,:) = simulate(tau_d(s,:), dt, N, ctrls{c}, A, B, E, C, ...
                            Kp, Ki, alpha_ff, gains{c}, Q_ref, w_ref, ...
                            R, Kt, Ke, Bd, Vmax, 1, L_obs);
    end
end

%% Metrics
% scenario 2: rise/OS/settle on the window before the disturbance fires,
% dip and recovery for the step itself. recovery band is 0.5% - the dips
% sit well inside the 2% settling band, which would read zero for everyone.
for s = 1:3
    if s == 2, i_end = find(t < t_dist(1), 1, 'last'); else, i_end = N; end

    fprintf('\n--- %s ---\n', names{s});
    fprintf('%-22s %8s %8s %10s %10s', 'ctrl','rise','OS[%]','settle','ss');
    if s == 2, fprintf(' %8s %10s','dip[%]','rec[0.5%]'); end
    fprintf('\n');

    for c = 1:nc
        y  = squeeze(Y(s,c,:))';
        yw = y(1:i_end);  tw = t(1:i_end);
        fprintf('%-22s %8.3f %8.2f %10.3f %10.2f', labels{c}, ...
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

%% Pole placement under load - the permanent offset, stated exactly
% mean over the last 0.5 s of the disturbance window, fully settled onto
% the shifted equilibrium the static prefilter cannot correct
y_pp = squeeze(Y(2,6,:))';
i_on = t >= t_dist(2)-0.5 & t < t_dist(2);
fprintf('\nPP under 4 mNm load: %.3f m^3/h (%.2f%% below target)\n', ...
        mean(y_pp(i_on)), (Q_target - mean(y_pp(i_on)))/Q_target*100);

%% dt convergence check - results should not move
y1 = simulate(tau_d(1,:), dt, N, 'LQR', A,B,E,C, Kp,Ki,alpha_ff,K_lqr, ...
              Q_ref,w_ref,R,Kt,Ke,Bd,Vmax);
dt2 = dt/4;  t2 = 0:dt2:T_sim;  N2 = numel(t2);
y2 = simulate(zeros(1,N2), dt2, N2, 'LQR', A,B,E,C, Kp,Ki,alpha_ff,K_lqr, ...
              Q_ref,w_ref,R,Kt,Ke,Bd,Vmax);
fprintf('\ndt check (LQR, nominal): rise %.4f @ dt=%.0e, %.4f @ dt=%.0e\n', ...
        riseTime(y1,t,Q_target), dt, riseTime(y2,t2,Q_target), dt2);

%% Feedforward sensitivity to load-measurement error
% the alpha sweep below scales the whole feedforward term - that is a
% feedforward-AUTHORITY sweep. this tests the measured-load assumption
% itself: the feedforward sees tau*(1+eps) while the plant sees tau.
fprintf('\nFF load-measurement error (filter clogging):\n%8s %8s %10s\n', ...
        'eps','dip[%]','rec[0.5%]');
eps_ff = [-0.2 0 0.2];
Yff = zeros(numel(eps_ff), N);
for i = 1:numel(eps_ff)
    Yff(i,:) = simulate(tau_d(2,:), dt, N, 'FF', A,B,E,C, Kp,Ki,alpha_ff,[], ...
                        Q_ref,w_ref,R,Kt,Ke,Bd,Vmax, 1+eps_ff(i));
    [d,rec] = dipMetrics(Yff(i,:),t,Q_target,t_dist,0.005);
    fprintf('%8.2f %8.3f %10.3f\n', eps_ff(i), d, rec);
end

%% PI gain grid - shows the baseline choice was deliberate
fprintf('\nPI gain grid (nominal load):\n%8s %9s %8s %8s %9s\n', ...
        'Kp','Ki','rise','OS[%]','settle');
for Kp_t = [2000 4000 8000]
    for ratio = [4 8 16]
        Ki_t = Kp_t*ratio;
        yg = simulate(tau_d(1,:), dt, N, 'PI', A,B,E,C, Kp_t,Ki_t, ...
                      alpha_ff,[], Q_ref,w_ref,R,Kt,Ke,Bd,Vmax);
        fprintf('%8d %9d %8.4f %8.3f %9.4f\n', Kp_t, Ki_t, ...
                riseTime(yg,t,Q_target), (max(yg)-Q_target)/Q_target*100, ...
                settleTime(yg,t,Q_target,0.02));
    end
end

%% Parameter sensitivity - J and Bd at +/-50%, controllers held fixed
% rise from the pre-disturbance window, dip from the +4 mNm step, so one
% run per corner backs both the rise-span and dip-ordering claims
fprintf('\nSensitivity corners (controllers at nominal design):\n');
fprintf('%6s %6s %9s %9s %9s %9s %9s %9s\n', 'J','Bd', ...
        'PI rise','FF rise','LQR rise','PI dip','FF dip','LQR dip');
for fJ = [0.5 1.5]
    for fB = [0.5 1.5]
        Jc = J*fJ;  Bc = Bd*fB;
        Ac = [-R/L, -Ke/L; Kt/Jc, -Bc/Jc];
        Ec = [0; -1/Jc];
        rr = zeros(1,3);  dd = zeros(1,3);
        modes = {'PI','FF','LQR'};
        for mi = 1:3
            yc = simulate(tau_d(2,:), dt, N, modes{mi}, Ac,B,Ec,C, ...
                          Kp,Ki,alpha_ff,K_lqr, Q_ref,w_ref, ...
                          R,Kt,Ke,Bc,Vmax);
            i_pre  = find(t < t_dist(1), 1, 'last');
            rr(mi) = riseTime(yc(1:i_pre), t(1:i_pre), Q_target);
            dd(mi) = dipMetrics(yc, t, Q_target, t_dist, 0.005);
        end
        fprintf('%6.2f %6.2f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n', ...
                fJ, fB, rr, dd);
    end
end

%% Observer pole placement sweep
% on the nominal scenario the sweep is degenerate by construction: the
% observer has the exact model, the same initial state, and the same
% delivered input u_sat, so the innovation is identically zero and L_obs
% never acts. pole choice only matters when the estimate is wrong, and
% the only error source here is the unmodeled load disturbance - so the
% sweep runs the filter-clogging step. all candidates respect the Euler
% limit (|p| < 2/dt = 1e5); the last pair is deliberately pathological:
% one pole slower than the controller's dominant pair, one 1.5x its fastest
fprintf('\nObserver pole placement (filter clogging):\n%22s %8s %10s\n', ...
        'poles','dip[%]','rec[0.5%]');
for pd = {[-300 -3000], [-800 -8000], [-2000 -20000], [-300 -1.5*40612]}
    Lo = place(A', C', pd{1})';
    yo = simulate(tau_d(2,:), dt, N, 'OBS', A,B,E,C, Kp,Ki,alpha_ff, ...
                  K_lqr, Q_ref,w_ref,R,Kt,Ke,Bd,Vmax, 1, Lo);
    [d,rec] = dipMetrics(yo,t,Q_target,t_dist,0.005);
    fprintf('%10.0f %11.0f %8.3f %10.3f\n', pd{1}, d, rec);
end

%% Plots
if ~exist('figs','dir'), mkdir('figs'); end
cols = {'b','g','r--','m-.','c','k--'};

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
% y-limits from the tracking controllers (all but the naive LQR, which
% sits near zero flow and just exits the frame - that is the point);
% including PP keeps its permanent post-step offset in frame
zwin = t >= 1.3 & t <= 2.3;
yz = Y(2,[1 2 3 5 6],zwin);
xlim([1.3 2.3]); ylim([min(yz(:))-0.3, max(yz(:))+0.3]);
ylabel('Flow [m^3/h]'); xlabel('Time [s]');
title('Load step at t = 1.5 s (zoom)');

delete(findall(f3,'Type','AxesToolbar'));
exportgraphics(f3, 'figs/zoom.png', 'Resolution', 150, 'BackgroundColor','white');

%% Feedforward measurement error vs the feedback-only LQR
% the point of this figure: the eps = -0.2 feedforward curve lands on top
% of the LQR's dip, so FF's advantage is exactly its load knowledge
f4 = figure('Color','w','Position',[50 50 700 420]);
if isprop(f4,'Theme'), f4.Theme = 'light'; end
hold on; grid on;
plot(t, Yff(2,:), 'g',   'LineWidth',1.4);
plot(t, Yff(1,:), 'g--', 'LineWidth',1.4);
plot(t, Yff(3,:), 'g:',  'LineWidth',1.4);
plot(t, squeeze(Y(2,3,:)), 'r--', 'LineWidth',1.4);
yline(Q_target,'k:','Target');
zw = t >= 1.4 & t <= 2.2;
yy = [Yff(:,zw); squeeze(Y(2,3,zw)).'];
xlim([1.4 2.2]); ylim([min(yy(:))-0.1, max(yy(:))+0.1]);
xlabel('Time [s]'); ylabel('Flow [m^3/h]');
title('Load step: feedforward seeing \tau(1+\epsilon) vs LQR');
legend({'FF, \epsilon = 0','FF, \epsilon = -0.2','FF, \epsilon = +0.2', ...
        'LQR (Bryson)'}, 'Location','SouthEast');
delete(findall(f4,'Type','AxesToolbar'));
exportgraphics(f4, 'figs/ff_error.png', 'Resolution', 150, 'BackgroundColor','white');

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
                      Q_ref, w_ref, R, Kt, Ke, Bd, Vmax, ff_scale, L_obs)
    if nargin < 20, ff_scale = 1; end   % load-measurement error, FF only
    if nargin < 21, L_obs = []; end     % observer gain, OBS only
    x = [0;0]; xhat = [0;0]; xi = 0; y = zeros(1,N);

    if strcmp(mode,'PP')
        % static prefilter inverting the nominal closed-loop DC gain;
        % the K_lqr slot carries the 1x2 state-feedback gain K_pp here
        G_pp = -1/(C/(A - B*K_lqr)*B);
    end

    for k = 1:N-1
        tau    = tau_d(k);
        y_meas = C*x;               % the only measurement available
        e      = Q_ref - y_meas;

        switch mode
            case 'PI'
                u = Kp*e + Ki*xi;
            case 'FF'
                % inverse dynamics at the reference; nominal drag Bd*w_ref
                % plus the measured disturbance torque
                V_ff = R*((Bd*w_ref + ff_scale*tau)/Kt) + Ke*w_ref;
                u = Kp*e + Ki*xi + alpha_ff*V_ff;
            case 'LQR'
                u = -K_lqr*[x; xi];
            case 'OBS'
                u = -K_lqr*[xhat; xi];      % control from the ESTIMATE
            case 'PP'
                u = -K_lqr*x + G_pp*Q_ref;  % no integral state at all
        end

        u_sat = min(Vmax, max(0,u));

        % conditional integration - without this the integrator charges
        % through the startup clamp and dumps it as overshoot on release
        if ~(u_sat ~= u && e*(u - u_sat) > 0)
            xi = xi + e*dt;                 % xi uses the real measurement
        end

        if strcmp(mode,'OBS')
            % driven by u_sat (what the actuator delivered, not what was
            % demanded) and with no E*tau term - the observer does not
            % know the disturbance, which is the honest setup
            innov = y_meas - C*xhat;
            xhat  = xhat + (A*xhat + B*u_sat + L_obs*innov)*dt;
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


function s = poleStr(p)
    s = '';
    for z = p.'
        if imag(z) == 0
            s = [s sprintf('%.6g, ', real(z))]; %#ok<AGROW>
        else
            s = [s sprintf('%.6g%+.6gj, ', real(z), imag(z))]; %#ok<AGROW>
        end
    end
    s = s(1:end-2);
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
