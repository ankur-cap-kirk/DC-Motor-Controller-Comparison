% Duct fan flow control (speed-referenced) - PI vs feedforward PI vs LQR
% with integral action, plus an observer-based LQR and a pole-placement
% strawman.
%
% Design: tangent-linearised model about the operating point, in deviation
% coordinates, discretised at the controller sample rate. Simulation: the
% nonlinear plant (w^2 aero load, Coulomb friction) under RK4, a 10 kHz
% discrete controller with input saturation, and sensor noise.
%
% Motor: FAULHABER 2342 S 024 CR, parameters from the series 2342...CR
% datasheet (edition 2026-03-23). Load: 80 mm axial impeller on the motor
% shaft. Flow is inferred from shaft speed by the fan affinity law at the
% nominal duct restriction - every controller here regulates speed-
% referenced flow, and restriction changes enter only as shaft-torque
% disturbances. See README, "What is actually being controlled".
%
% A. Guruprasad, June 2026. Revised November 2026 - see README, Revision notes.

clear; clc; close all;

%% Motor - Faulhaber 2342S024CR, 24 V
% datasheet: R 7.14 ohm, L 264 uH, kM 28.6 mNm/A, kn 334 rpm/V,
% rotor inertia 5.8 gcm^2, rated 20.5 mNm / 0.846 A / 5660 rpm
R  = 7.14;          L  = 264e-6;
Kt = 28.6e-3;       Ke = 60/(2*pi*334);     % 334 rpm/V -> V.s/rad
J_rot = 5.8e-7;
Vmax  = 24;         I_rated = 0.846;        I_stall = Vmax/R;

% Kt and Ke are the same constant in SI - if they disagree, one of the
% conversions above is wrong
assert(abs(Kt - Ke)/Kt < 0.01, 'Kt/Ke mismatch - unit conversion error');
fprintf('Kt = %.4f N.m/A vs Ke = %.4f V.s/rad  (agree to %.2f%%)\n', ...
        Kt, Ke, abs(Kt - Ke)/Kt*100);

%% Load - 80 mm axial impeller, direct drive. Nonlinear static map.
% impeller inertia from mass and radius: 20 g at 40 mm tip radius as a
% solid disc, J = 0.5*m*r^2. estimate, and a high one - real impellers
% carry most of their mass near the hub
m_imp = 0.020;  r_imp = 0.040;
J = J_rot + 0.5*m_imp*r_imp^2;

% operating point: 4000 rpm pushing 85 m^3/h through the nominal duct
w_op  = 4000*2*pi/60;
Q_op  = 85/3600;
K_out = Q_op/w_op;              % [m^3/s per rad/s]

% aero torque follows the fan law, tau = k_a*w^2, pinned to 8 mNm at w_op
% (~3.4 W shaft power, sane for an 80 mm fan at 4000 rpm). the 0.99 mNm
% datasheet friction torque is Coulomb, not viscous. estimates, not data.
k_a   = 8.0e-3/w_op^2;          % [N.m/(rad/s)^2]
tau_c = 0.99e-3;                % [N.m]

% stay inside the continuous ratings at the operating point
tau_ss = k_a*w_op^2 + tau_c;  I_ss = tau_ss/Kt;  V_ss = R*I_ss + Ke*w_op;
x_op   = [I_ss; w_op];
fprintf('operating point: %.1f mNm, %.3f A, %.2f V  (rated: 20.5 mNm, %.3f A, %d V)\n', ...
        tau_ss*1e3, I_ss, V_ss, I_rated, Vmax);
assert(tau_ss < 20.5e-3 && I_ss < I_rated && V_ss < Vmax, 'operating point exceeds ratings');

Q_target = 85;                  % [m^3/h] for display

%% Linear design model - TANGENT linearisation about the operating point
% x = [i; w], u = terminal voltage, tau_d = load torque deviation [N.m],
% y = flow [m^3/s]. the small-signal aero damping is the slope of k_a*w^2
% at w_op, i.e. 2*k_a*w_op - twice the secant tau/w. Coulomb friction
% contributes no slope. everything model-based below works in deviations
% from (x_op, V_ss).
B_lin = 2*k_a*w_op;
A = [-R/L, -Ke/L; Kt/J, -B_lin/J];
B = [1/L; 0];
E = [0; -1/J];
C = [0, K_out];

fprintf('open-loop poles: %s\n', poleStr(eig(A)));
fprintf('rank(ctrb) = %d, rank(obsv) = %d\n', rank(ctrb(A,B)), rank(obsv(A,C)));
fprintf('aero damping: tangent %.3e N.m.s (secant tau/w would be %.3e, %.0f%% low)\n', ...
        B_lin, tau_ss/w_op, (1 - (tau_ss/w_op)/B_lin)*100);

%% Controller sample time and plant integration step
% the electrical pole (~27 krad/s) is far above any sensible controller
% rate, so the controllers are designed in discrete time at Ts on the ZOH
% model rather than pretending to be continuous. the plant itself is
% integrated with RK4 at a finer substep, input held over each Ts.
Ts   = 1e-4;  nsub = 5;  dt = Ts/nsub;
sysd = c2d(ss(A,B,C,0), Ts, 'zoh');  Ad = sysd.A;  Bd = sysd.B;
assert(max(abs(eig(A)))*dt < 2.5, 'RK4 substep too coarse for the electrical pole');
fprintf('controller Ts = %.0e s (%d kHz); plant RK4 substep %.0e s\n', Ts, round(1e-3/Ts), dt);

%% LQR with integral state - deviation coordinates, discrete design
% xi = sum Ts*(Q_cmd - y). the regulator drives DEVIATIONS to zero, so the
% equilibrium is the origin and the integrator only has to hold the
% operating-point voltage V_ss (plus any disturbance bias) - not fight a
% speed-feedback term demanding w -> 0.
Aa = [Ad, zeros(2,1); -Ts*C, 1];
Ba = [Bd; 0];
assert(rank(ctrb(Aa,Ba)) == 3, 'augmented pair not controllable - check C');

% Bryson weighting. i ~ O(1), w ~ O(400), xi ~ O(1e-5): Q = eye(3) lets
% speed dominate and leaves the integral state effectively unweighted
% (kept below as K_naive to show exactly that). Maxima:
%   max_i   3 A     transient excursion, stall is 3.36 A
%   max_w   40      ~10% of operating speed
%   max_xi  2.8e-5  ~1 m^3/h steady error corrected in ~0.1 s
%   max_u   24      the rail
max_i = 3;  max_w = 40;  max_xi = 2.8e-5;  max_u = 24;
Qw = diag([1/max_i^2, 1/max_w^2, 1/max_xi^2]);
Rw = 1/max_u^2;

K_lqr = dlqr(Aa, Ba, Qw, Rw);
p_lqr = log(eig(Aa - Ba*K_lqr))/Ts;           % s-plane equivalents
fprintf('K_lqr = [%.4f %.4f %.4e]\n', K_lqr);
fprintf('LQR closed-loop poles (s-equiv): %s\n', poleStr(p_lqr));
fprintf('integrator holds the operating point at xi = %.2e  (Bryson max_xi = %.1e)\n', ...
        V_ss/(-K_lqr(3)), max_xi);

% the naive design, kept in the comparison on purpose
K_naive = dlqr(Aa, Ba, eye(3), 1);
p_naive = log(eig(Aa - Ba*K_naive))/Ts;
fprintf('naive LQR closed-loop poles (s-equiv): %s\n', poleStr(p_naive));

%% Luenberger observer - estimates current from the flow measurement
% (A,C) is observable, so full state feedback does not require a current
% sensor. discrete observer in deviation coordinates, driven by the
% DELIVERED input and the flow measurement only - it is given no knowledge
% of the load disturbance. the electrical mode is left where the plant
% puts it (moving a 27 krad/s pole through a speed measurement costs
% enormous gain and buys nothing); only the mechanical-estimate pole is a
% design choice, and it is a speed/noise tradeoff - see the sweep below.
p_el  = min(real(eig(A)));
p_obs = [-2000, p_el];
L_obs = place(Ad', C', exp(p_obs*Ts))';
fprintf('observer poles (s-equiv): %s ; L_obs = [%.3e %.3e]\n', ...
        poleStr(log(eig(Ad - L_obs*C))/Ts), L_obs);

%% Pole placement with static prefilter - the no-integrator strawman
% poles taken from the LQR's own closed loop (fastest real part, dominant
% real part), u = V_ss - K_pp*dx + G_pp*(Q_cmd - Q_op) with G_pp inverting
% the nominal DC gain. exact on the nominal plant; no mechanism to notice
% the plant changed, so any load deviation becomes a permanent offset.
pr   = sort(real(p_lqr));
p_pp = [pr(1), pr(3)];
K_pp = place(Ad, Bd, exp(p_pp*Ts));
G_pp = 1/(C/(eye(2) - Ad + Bd*K_pp)*Bd);
fprintf('K_pp = [%.4f %.4f] at poles %s ; G_pp = %.4e\n', K_pp, poleStr(p_pp), G_pp);

%% PI design - analytic, from the reduced first-order plant
% Timescale separation (electrical pole ~3000x faster than mechanical)
% justifies di/dt -> 0, giving P(s) = kdc/(tau*s + 1) with the back-EMF
% damping Kt*Ke/R folded into the effective friction. PI then yields an
% exactly second-order closed loop:
%     Kp = (2*zeta*wn*tau - 1)/kdc
%     Ki = wn^2*tau/kdc
% The PI zero at -Ki/Kp lands next to the slow closed-loop pole, so the
% linear response is effectively first order at the fast pole. This
% ignores input saturation, which binds hard during startup - the small-
% step scenario is where the design actually gets tested.
B_eff = B_lin + Kt*Ke/R;
tau_m = J/B_eff;
kdc   = K_out*Kt/(R*B_eff);

zeta_d = 1.50;  wn_d = 20.9;   % deliberately overdamped baseline
Kp = (2*zeta_d*wn_d*tau_m - 1)/kdc;
Ki = wn_d^2*tau_m/kdc;
p_pi = roots([1, 2*zeta_d*wn_d, wn_d^2]);
fprintf('PI design: tau_m = %.4f s, kdc = %.3e; zeta = %.2f, wn = %.1f -> Kp = %.0f V/(m^3/s), Ki = %.0f V/m^3\n', ...
        tau_m, kdc, zeta_d, wn_d, Kp, Ki);
fprintf('   closed-loop poles %s, PI zero at %.2f (near-cancellation of the slow pole)\n', ...
        poleStr(p_pi), -Ki/Kp);

alpha_ff = 0.86;    % knee of the sweep below: rise time is on the
                    % saturation floor there, overshoot starts ~0.90

%% Loop margins - the discrete loops as implemented
Cpi = tf([Kp, Ki*Ts - Kp], [1, -1], Ts);      % Kp + Ki*Ts/(z-1)
[gm, pm, ~, wc] = margin(Cpi*sysd);
fprintf('PI loop:  GM %.1f dB, PM %.1f deg, crossover %.1f rad/s\n', 20*log10(gm), pm, wc);
[gm, pm, ~, wc] = margin(ss(Aa, Ba, K_lqr, 0, Ts));
fprintf('LQR loop (broken at plant input): GM %.1f dB, PM %.1f deg, crossover %.1f rad/s\n', ...
        20*log10(gm), pm, wc);

%% Sensor noise - white at the sample rate
sig_y = 0.25/3600;   % flow sensor, 1-sigma [m^3/s]: 0.25 m^3/h, 0.3% of setpoint
sig_i = 0.010;       % current sensor, 1-sigma [A]: full-state designs only
fprintf('sensor noise: flow %.2f m^3/h, current %.0f mA (1-sigma, white at Ts)\n', sig_y*3600, sig_i*1e3);

%% Plant (truth) and model (what the controllers know)
P = struct('R',R,'L',L,'Kt',Kt,'Ke',Ke,'J',J,'k_a',k_a,'tau_c',tau_c,'Vmax',Vmax);
M = struct('Ad',Ad,'Bd',Bd,'C',C,'x_op',x_op,'V_ss',V_ss,'Q_op',Q_op,'K_out',K_out, ...
           'R',R,'Kt',Kt,'Ke',Ke,'k_a',k_a,'tau_c',tau_c,'Ts',Ts,'nsub',nsub, ...
           'sig_y',sig_y,'sig_i',sig_i);

ctl = { mkCtrl('PI',  'Standard PI',            'Kp',Kp,'Ki',Ki), ...
        mkCtrl('FF',  'Feedforward PI',         'Kp',Kp,'Ki',Ki,'alpha',alpha_ff), ...
        mkCtrl('LQR', 'LQR (Bryson)',           'K',K_lqr), ...
        mkCtrl('LQR', 'LQR (Q = I)',            'K',K_naive), ...
        mkCtrl('OBS', 'LQR + observer',         'K',K_lqr,'L',L_obs), ...
        mkCtrl('PP',  'Pole place + prefilter', 'K',K_pp,'G',G_pp) };
nc     = numel(ctl);
labels = cellfun(@(c) c.label, ctl, 'UniformOutput', false);

%% Scenarios
% 1. small step 80 -> 85 m^3/h from the settled 80 state: no saturation,
%    this is where the linear designs are actually compared
% 2. startup from rest to 85: a saturation study - every controller sits
%    on the 24 V rail at stall current, and differs only in how it leaves
% 3. filter clogging: +4 mNm shaft-torque step at t = 1.5 s, held to 2.75 s
% 4. gusty inlet: continuously varying torque - persistent excitation
Q_small = 80/3600;
sc = cell(1,4);

sc{1} = mkScen('Small step 80 -> 85', 1.0, Ts, Q_op);
sc{1}.Q_cmd(:) = Q_small;  sc{1}.Q_cmd(sc{1}.t >= 0.2) = Q_op;
sc{1}.x0 = equilibrium(Q_small, M);  sc{1}.preload = true;  sc{1}.t_step = 0.2;

sc{2} = mkScen('Startup (nominal load)', 1.5, Ts, Q_op);

sc{3} = mkScen('Filter clogging', 3.0, Ts, Q_op);
t_dist = [1.5 2.75];
sc{3}.tau_d(sc{3}.t > t_dist(1) & sc{3}.t < t_dist(2)) = 4e-3;

sc{4} = mkScen('Gusty inlet', 4.0, Ts, Q_op);
rng(1);
gust_lp = Ts/0.05;      % 50 ms first-order filter on the gust noise
g = filter(gust_lp, [1, -(1-gust_lp)], randn(1, sc{4}.N));
sc{4}.tau_d = 1.5e-3*sin(3*sc{4}.t) + 6e-4*(g/std(g));

ns = numel(sc);

%% Run
Y = cell(ns, nc);
for s = 1:ns
    for c = 1:nc
        Y{s,c} = simulate(P, M, ctl{c}, sc{s});
    end
end

%% Metrics - computed on the TRUE flow; the controllers only ever see the
% noisy measurement. u-std is the actuator jitter in steady state, the cost
% of a controller's noise sensitivity.

% 1. small step: rise/OS/settle relative to the 5 m^3/h step
fprintf('\n--- %s ---  (rise 10-90%% of step, OS %% of step, settle to 5%% of step)\n', sc{1}.name);
fprintf('%-22s %8s %8s %10s %9s %8s\n', 'ctrl','rise','OS[%]','settle','peak u','u-std');
o1 = Y{1,1}; k0 = find(o1.t >= sc{1}.t_step, 1);
for c = 1:nc
    o = Y{1,c};  ys = o.y(k0:end);  ts = o.t(k0:end) - o.t(k0);
    y0 = Q_small*3600;  dy = Q_target - y0;
    fprintf('%-22s %8.4f %8s %10.4f %9.1f %8.3f\n', labels{c}, ...
            riseTime(ys, ts, y0, Q_target), osStr(ys, y0, Q_target, dy), ...
            settleTime(ys, ts, Q_target, 0.05*dy), max(o.u), std(o.u(o.t > 0.7)));
end

% 2. startup: the saturation study - report what the rail did
fprintf('\n--- %s ---  (rise 10-90%%, settle to 2%%, SS = mean of last 0.5 s)\n', sc{2}.name);
fprintf('%-22s %8s %8s %10s %8s %8s %8s %8s\n', 'ctrl','rise','OS[%]','settle','ss','peak i','rail[s]','u-std');
for c = 1:nc
    o = Y{2,c};
    rail = o.t(find(o.u >= Vmax - 1e-9, 1, 'last'));  if isempty(rail), rail = 0; end
    fprintf('%-22s %8.3f %8s %10.3f %8.2f %8.2f %8.3f %8.3f\n', labels{c}, ...
            riseTime(o.y, o.t, 0, Q_target), osStr(o.y, 0, Q_target, Q_target), ...
            settleTime(o.y, o.t, Q_target, 0.02*Q_target), mean(o.y(o.t > o.t(end) - 0.5)), ...
            max(o.i), rail, std(o.u(o.t > o.t(end) - 0.5)));
end

% 3. load step: dip and recovery to 0.5% (the dips sit inside the 2%
% settling band, which would read zero for everyone); SS measured WHILE
% the load is on, so a permanent offset shows up where it belongs
fprintf('\n--- %s ---  (dip and recovery to 0.5%% for the +4 mNm step, SS = mean over t in [2.25, 2.75] s)\n', sc{3}.name);
fprintf('%-22s %8s %10s %12s\n', 'ctrl','dip[%]','rec[0.5%]','ss under load');
i_load = sc{3}.t >= t_dist(2) - 0.5 & sc{3}.t < t_dist(2);
for c = 1:nc
    o = Y{3,c};
    pre = o.y(find(o.t < t_dist(1), 1, 'last'));
    if abs(pre - Q_target) < 0.02*Q_target
        [d, rec] = dipMetrics(o.y, o.t, Q_target, t_dist, 0.005);
        fprintf('%-22s %8.2f %10.3f %12.3f\n', labels{c}, d, rec, mean(o.y(i_load)));
    else
        fprintf('%-22s %8s %10s %12.3f   (never settled before the step)\n', labels{c}, '-', '-', mean(o.y(i_load)));
    end
end

% 4. gust: tracking error statistics once the startup transient is over
fprintf('\n--- %s ---  (t > 1 s)\n', sc{4}.name);
fprintf('%-22s %10s %10s %8s\n', 'ctrl','RMS err','max |err|','u-std');
for c = 1:nc
    o = Y{4,c};  m = o.t > 1;  err = o.y(m) - Q_target;
    fprintf('%-22s %10.3f %10.3f %8.3f\n', labels{c}, rms(err), max(abs(err)), std(o.u(m)));
end

%% Integration check - results should not move with the substep
M2 = M;  M2.nsub = 20;
oa = Y{2,3};  ob = simulate(P, M2, ctl{3}, sc{2});
fprintf('\nsubstep check (LQR, startup): rise %.4f @ dt=%.0e, %.4f @ dt=%.0e\n', ...
        riseTime(oa.y, oa.t, 0, Q_target), Ts/M.nsub, riseTime(ob.y, ob.t, 0, Q_target), Ts/M2.nsub);

%% Feedforward sensitivity to load-measurement error
% the alpha sweep below scales the whole feedforward term - that is a
% feedforward-AUTHORITY sweep. this tests the measured-load assumption
% itself: the feedforward sees tau*(1+eps) while the plant sees tau.
fprintf('\nFF load-measurement error (filter clogging):\n%8s %8s %10s\n', 'eps','dip[%]','rec[0.5%]');
eps_ff = [-0.2 0 0.2];
Yff = cell(1, numel(eps_ff));
for i = 1:numel(eps_ff)
    cf = ctl{2};  cf.ffs = 1 + eps_ff(i);
    Yff{i} = simulate(P, M, cf, sc{3});
    [d, rec] = dipMetrics(Yff{i}.y, Yff{i}.t, Q_target, t_dist, 0.005);
    fprintf('%8.2f %8.3f %10.3f\n', eps_ff(i), d, rec);
end

%% PI gain grid - small step for the linear design, load step for rejection
fprintf('\nPI gain grid (small step 80 -> 85 | filter clogging):\n%8s %9s %8s %8s %9s %8s | %8s %10s\n', ...
        'Kp','Ki','rise','OS[%]','settle','u-std','dip[%]','rec[0.5%]');
for Kp_t = [2000 4000 8000]
    for ratio = [4 8 16]
        Ki_t = Kp_t*ratio;
        cg = mkCtrl('PI', 'grid', 'Kp', Kp_t, 'Ki', Ki_t);
        o  = simulate(P, M, cg, sc{1});
        ys = o.y(k0:end);  ts = o.t(k0:end) - o.t(k0);
        o3 = simulate(P, M, cg, sc{3});
        [d, rec] = dipMetrics(o3.y, o3.t, Q_target, t_dist, 0.005);
        fprintf('%8d %9d %8.4f %8s %9.4f %8.3f | %8.2f %10.3f\n', Kp_t, Ki_t, ...
                riseTime(ys, ts, y0, Q_target), osStr(ys, y0, Q_target, dy), ...
                settleTime(ys, ts, Q_target, 0.05*dy), std(o.u(o.t > 0.7)), d, rec);
    end
end

%% Parameter sensitivity - J and k_a at +/-50%, controllers held fixed
% the plant is perturbed; every controller (and the observer's model)
% keeps the nominal design. rise from the pre-disturbance window, dip from
% the +4 mNm step. the observer is the one design with a model to be
% wrong, so it is in the table.
fprintf('\nSensitivity corners (controllers at nominal design, load step):\n');
fprintf('%6s %6s | %7s %7s %7s %7s | %7s %7s %7s %7s\n', 'J','k_a', ...
        'PI','FF','LQR','OBS', 'PI','FF','LQR','OBS');
fprintf('%13s | %31s | %31s\n', '', 'rise [s]', 'dip [%]');
cidx = [1 2 3 5];
for fJ = [0.5 1.5]
    for fK = [0.5 1.5]
        Pc = P;  Pc.J = J*fJ;  Pc.k_a = k_a*fK;
        rr = zeros(1,4);  dd = zeros(1,4);
        for m = 1:4
            o = simulate(Pc, M, ctl{cidx(m)}, sc{3});
            ipre  = o.t < t_dist(1);
            rr(m) = riseTime(o.y(ipre), o.t(ipre), 0, Q_target);
            dd(m) = dipMetrics(o.y, o.t, Q_target, t_dist, 0.005);
        end
        fprintf('%6.2f %6.2f | %7.3f %7.3f %7.3f %7.3f | %7.3f %7.3f %7.3f %7.3f\n', fJ, fK, rr, dd);
    end
end

%% Observer pole placement sweep - speed against noise
% on an undisturbed run with an exact model the innovation is zero and the
% poles are invisible, so the sweep uses the load step. the electrical
% pole stays at the plant's; the mechanical-estimate pole is swept. faster
% tracks the disturbance sooner and passes more sensor noise to the
% actuator - u-std is that cost. the last entry slows the electrical mode
% as well, to show what that costs.
fprintf('\nObserver pole placement (filter clogging):\n%22s %8s %10s %8s\n', ...
        'poles (s-equiv)','dip[%]','rec[0.5%]','u-std');
for pd = {[-500 p_el], [-2000 p_el], [-8000 p_el], [-20000 p_el], [-2000 -20000]}
    Lo = place(Ad', C', exp(pd{1}*Ts))';
    co = mkCtrl('OBS', 'sweep', 'K', K_lqr, 'L', Lo);
    o  = simulate(P, M, co, sc{3});
    [d, rec] = dipMetrics(o.y, o.t, Q_target, t_dist, 0.005);
    fprintf('%10.0f %11.0f %8.3f %10.3f %8.3f\n', pd{1}, d, rec, std(o.u(o.t > 1 & o.t < 1.5)));
end

%% Plots
if ~exist('figs','dir'), mkdir('figs'); end
cols = {'b','g','r--','m-.','c','k--'};

% all scenarios (startup, load, gust) with their disturbance traces
f1 = figure('Color','w','Position',[50 50 1100 850]);
if isprop(f1,'Theme'), f1.Theme = 'light'; end     % R2025b batch defaults to dark
for s = 2:4
    r = s - 1;
    subplot(3,2,2*r-1); hold on; grid on;
    for c = 1:nc, plot(Y{s,c}.t, Y{s,c}.y, cols{c}, 'LineWidth',1.4); end
    yline(Q_target,'k:','Target');
    ylabel('Flow [m^3/h]'); ylim([0 110]); title(sc{s}.name);
    if r == 1, legend(labels,'Location','SouthEast'); end
    if r == 3, xlabel('Time [s]'); end

    subplot(3,2,2*r); grid on;
    area(sc{s}.t, sc{s}.tau_d*1e3, 'FaceColor', [.85 .85 .85]);
    ylabel('\tau_d [mNm]'); ylim([-3 6]);
    if r == 3, xlabel('Time [s]'); end
end
delete(findall(f1,'Type','AxesToolbar'));
exportgraphics(f1, 'figs/comparison.png', 'Resolution', 150, 'BackgroundColor','white');

% zoomed views - startup and the load-step dip
f3 = figure('Color','w','Position',[50 50 1100 420]);
if isprop(f3,'Theme'), f3.Theme = 'light'; end
subplot(1,2,1); hold on; grid on;
for c = 1:nc, plot(Y{2,c}.t, Y{2,c}.y, cols{c}, 'LineWidth',1.4); end
yline(Q_target,'k:','Target');
xlim([0 0.6]); ylabel('Flow [m^3/h]'); xlabel('Time [s]');
title('Startup (nominal load)'); legend(labels,'Location','SouthEast');

subplot(1,2,2); hold on; grid on;
for c = 1:nc, plot(Y{3,c}.t, Y{3,c}.y, cols{c}, 'LineWidth',1.4); end
yline(Q_target,'k:','Target');
% y-limits from the designs that hold the setpoint; the naive LQR sits a
% few percent low and exits the frame, which is the point
zwin = sc{3}.t >= 1.3 & sc{3}.t <= 2.3;
yz = cellfun(@(o) [min(o.y(zwin)), max(o.y(zwin))], Y(3,[1 2 3 5 6]), 'UniformOutput', false);
yz = vertcat(yz{:});
xlim([1.3 2.3]); ylim([min(yz(:,1))-0.3, max(yz(:,2))+0.3]);
ylabel('Flow [m^3/h]'); xlabel('Time [s]'); title('Load step at t = 1.5 s (zoom)');
delete(findall(f3,'Type','AxesToolbar'));
exportgraphics(f3, 'figs/zoom.png', 'Resolution', 150, 'BackgroundColor','white');

% small-signal step: the linear comparison, with the control effort
f5 = figure('Color','w','Position',[50 50 1100 620]);
if isprop(f5,'Theme'), f5.Theme = 'light'; end
subplot(2,1,1); hold on; grid on;
for c = 1:nc, plot(Y{1,c}.t, Y{1,c}.y, cols{c}, 'LineWidth',1.4); end
plot(Y{1,1}.t, Y{1,1}.Q_cmd, 'k:', 'LineWidth',1);
xlim([0.15 0.6]); ylim([79.5 86]); ylabel('Flow [m^3/h]');
title('Small step 80 \rightarrow 85 m^3/h from the settled state (no saturation)');
legend([labels, {'Command'}],'Location','SouthEast');
subplot(2,1,2); hold on; grid on;
for c = 1:nc, plot(Y{1,c}.t, Y{1,c}.u, cols{c}, 'LineWidth',1.2); end
yline(Vmax,'k:','24 V rail');
xlim([0.15 0.6]); ylabel('u [V]'); xlabel('Time [s]');
delete(findall(f5,'Type','AxesToolbar'));
exportgraphics(f5, 'figs/smallstep.png', 'Resolution', 150, 'BackgroundColor','white');

% control effort during startup: the rail and the current it draws
f6 = figure('Color','w','Position',[50 50 1100 620]);
if isprop(f6,'Theme'), f6.Theme = 'light'; end
subplot(2,1,1); hold on; grid on;
for c = 1:nc, plot(Y{2,c}.t, Y{2,c}.u, cols{c}, 'LineWidth',1.2); end
yline(V_ss,'k:','V_{ss}');
xlim([0 0.3]); ylim([0 26]); ylabel('u [V]'); title('Startup: terminal voltage and armature current');
legend(labels,'Location','NorthEast');
subplot(2,1,2); hold on; grid on;
for c = 1:nc, plot(Y{2,c}.t, Y{2,c}.i, cols{c}, 'LineWidth',1.2); end
yline(I_rated,'k:','rated'); yline(I_stall,'k--','stall');
xlim([0 0.3]); ylim([0 3.6]); ylabel('i [A]'); xlabel('Time [s]');
delete(findall(f6,'Type','AxesToolbar'));
exportgraphics(f6, 'figs/effort.png', 'Resolution', 150, 'BackgroundColor','white');

% feedforward measurement error vs the feedback-only LQR
f4 = figure('Color','w','Position',[50 50 700 420]);
if isprop(f4,'Theme'), f4.Theme = 'light'; end
hold on; grid on;
plot(Yff{2}.t, Yff{2}.y, 'g',   'LineWidth',1.4);
plot(Yff{1}.t, Yff{1}.y, 'g--', 'LineWidth',1.4);
plot(Yff{3}.t, Yff{3}.y, 'g:',  'LineWidth',1.4);
plot(Y{3,3}.t, Y{3,3}.y, 'r--', 'LineWidth',1.4);
yline(Q_target,'k:','Target');
zw = sc{3}.t >= 1.4 & sc{3}.t <= 2.2;
yy = [Yff{1}.y(zw); Yff{2}.y(zw); Yff{3}.y(zw); Y{3,3}.y(zw)];
xlim([1.4 2.2]); ylim([min(yy(:))-0.1, max(yy(:))+0.1]);
xlabel('Time [s]'); ylabel('Flow [m^3/h]');
title('Load step: feedforward seeing \tau(1+\epsilon) vs LQR');
legend({'FF, \epsilon = 0','FF, \epsilon = -0.2','FF, \epsilon = +0.2','LQR (Bryson)'}, 'Location','SouthEast');
delete(findall(f4,'Type','AxesToolbar'));
exportgraphics(f4, 'figs/ff_error.png', 'Resolution', 150, 'BackgroundColor','white');

%% Feedforward scaling sweep (startup)
a  = 0.5:0.05:1.0;
rt = zeros(size(a));  os = zeros(size(a));
for k = 1:numel(a)
    ca = ctl{2};  ca.alpha = a(k);
    o  = simulate(P, M, ca, sc{2});
    rt(k) = riseTime(o.y, o.t, 0, Q_target);
    os(k) = (max(o.y) - Q_target)/Q_target*100;
end
fprintf('\nfeedforward sweep (startup):\n%8s %8s %8s\n', 'alpha','rise','OS[%]');
fprintf('%8.2f %8.3f %8.2f\n', [a; rt; os]);

f2 = figure('Color','w');
if isprop(f2,'Theme'), f2.Theme = 'light'; end
yyaxis left;  plot(a, rt, 'LineWidth',1.5); ylabel('Rise time [s]');
yyaxis right; plot(a, os, 'LineWidth',1.5); ylabel('Overshoot [%]');
xlabel('\alpha'); grid on; title('Feedforward scaling tradeoff');
delete(findall(f2,'Type','AxesToolbar'));
exportgraphics(f2, 'figs/ff_sweep.png', 'Resolution', 150, 'BackgroundColor','white');


%% ---------------------------------------------------------------------
function c = mkCtrl(mode, label, varargin)
    c = struct('mode',mode, 'label',label, 'Kp',0, 'Ki',0, 'alpha',0, ...
               'ffs',1, 'K',[], 'L',[], 'G',0);
    for k = 1:2:numel(varargin), c.(varargin{k}) = varargin{k+1}; end
end


function s = mkScen(name, T, Ts, Q_op)
    s.name  = name;  s.t = 0:Ts:T;  s.N = numel(s.t);
    s.Q_cmd = Q_op*ones(1, s.N);  s.tau_d = zeros(1, s.N);
    s.x0 = [0; 0];  s.preload = false;  s.seed = 0;  s.t_step = 0;
end


function [x0, u0] = equilibrium(Q, M)
    % steady state of the NONLINEAR plant at flow Q, zero disturbance
    w0 = Q/M.K_out;
    i0 = (M.k_a*w0^2 + M.tau_c)/M.Kt;
    u0 = M.R*i0 + M.Ke*w0;
    x0 = [i0; w0];
end


function V = ffVolt(Qc, tau, ffs, M)
    % inverse dynamics at the commanded flow: nominal static load plus
    % the (scaled) measured disturbance torque
    w_ref = Qc/M.K_out;
    V = M.R*((M.k_a*w_ref^2 + M.tau_c + ffs*tau)/M.Kt) + M.Ke*w_ref;
end


function dx = plantDeriv(x, u, tau_d, R, L, Kt, Ke, J, k_a, tau_c)
    % the truth: w^2 aero load and Coulomb friction (smoothed at |w| < 0.5)
    i = x(1);  w = x(2);
    dx = [(u - R*i - Ke*w)/L;
          (Kt*i - k_a*w*abs(w) - tau_c*tanh(w/0.5) - tau_d)/J];
end


function o = simulate(P, M, c, s)
    R = P.R; L = P.L; Kt = P.Kt; Ke = P.Ke; J = P.J; k_a = P.k_a; tau_c = P.tau_c; Vmax = P.Vmax;
    Ad = M.Ad; Bd = M.Bd; C = M.C; x_op = M.x_op; V_ss = M.V_ss; Q_op = M.Q_op; K_out = M.K_out;
    Ts = M.Ts; nsub = M.nsub; dt = Ts/nsub;
    N = s.N;  x = s.x0;  dxh = [0; 0];  xi = 0;

    rng(s.seed);                        % same noise sequence for every controller
    ny = M.sig_y*randn(1, N);
    ni = M.sig_i*randn(1, N);

    % integrator preload: start holding the initial equilibrium exactly,
    % so a reference step tests the design, not the controller switch-on
    if s.preload
        [~, u0] = equilibrium(C*x, M);
        dx0 = x - x_op;
        switch c.mode
            case 'PI',  xi = u0/c.Ki;
            case 'FF',  xi = (u0 - c.alpha*ffVolt(s.Q_cmd(1), 0, c.ffs, M))/c.Ki;
            case {'LQR','OBS'}
                xi = (u0 + c.K(1:2)*dx0)/(-c.K(3));  dxh = dx0;
        end
    end

    y = zeros(1,N); u = zeros(1,N); iA = zeros(1,N); w = zeros(1,N); xh = zeros(2,N);
    for k = 1:N
        y(k) = C*x;  iA(k) = x(1);  w(k) = x(2);  xh(:,k) = dxh + x_op;

        ym = y(k) + ny(k);              % the flow measurement, noisy
        Qc = s.Q_cmd(k);  e = Qc - ym;

        switch c.mode
            case 'PI'
                uc = c.Kp*e + c.Ki*xi;
            case 'FF'
                uc = c.Kp*e + c.Ki*xi + c.alpha*ffVolt(Qc, s.tau_d(k), c.ffs, M);
            case 'LQR'                  % two sensors: current and speed
                dxm = [x(1) + ni(k) - x_op(1);  ym/K_out - x_op(2)];
                uc  = -c.K*[dxm; xi];
            case 'OBS'                  % one sensor: control from the ESTIMATE
                uc  = -c.K*[dxh; xi];
            case 'PP'                   % two sensors, no integral state
                dxm = [x(1) + ni(k) - x_op(1);  ym/K_out - x_op(2)];
                uc  = V_ss - c.K*dxm + c.G*(Qc - Q_op);
        end

        us = min(Vmax, max(0, uc));  u(k) = us;

        % conditional integration - without this the integrator charges
        % through the startup clamp and dumps it as overshoot on release.
        % the integrator is left at zero when the rail releases, which is
        % why startup rise times are not the linear designs' rise times.
        if ~(us ~= uc && e*(uc - us) > 0)
            xi = xi + Ts*e;
        end

        % observer: driven by the DELIVERED input and the measurement only;
        % no disturbance term - it does not know the load
        if strcmp(c.mode, 'OBS')
            dxh = Ad*dxh + Bd*(us - V_ss) + c.L*((ym - Q_op) - C*dxh);
        end

        % plant: RK4 over the control period with the input held
        td = s.tau_d(k);
        for m = 1:nsub
            a1 = plantDeriv(x,             us, td, R, L, Kt, Ke, J, k_a, tau_c);
            a2 = plantDeriv(x + 0.5*dt*a1, us, td, R, L, Kt, Ke, J, k_a, tau_c);
            a3 = plantDeriv(x + 0.5*dt*a2, us, td, R, L, Kt, Ke, J, k_a, tau_c);
            a4 = plantDeriv(x + dt*a3,     us, td, R, L, Kt, Ke, J, k_a, tau_c);
            x  = x + dt/6*(a1 + 2*a2 + 2*a3 + a4);
        end
    end

    o.t = s.t;  o.y = y*3600;  o.u = u;  o.i = iA;  o.w = w;  o.xhat = xh;
    o.Q_cmd = s.Q_cmd*3600;     % m^3/s -> m^3/h for display
end


function tr = riseTime(y, t, y0, y1)
    i10 = find(y >= y0 + 0.10*(y1 - y0), 1);
    i90 = find(y >= y0 + 0.90*(y1 - y0), 1);
    if isempty(i10) || isempty(i90), tr = NaN; else, tr = t(i90) - t(i10); end
end


function ts = settleTime(y, t, target, band)
    out = find(abs(y - target) > band);
    if isempty(out), ts = 0; else, ts = t(out(end)); end
end


function s = osStr(y, y0, y1, denom)
    % overshoot past y1 as a percentage of denom; n/a for a trace that
    % never gets within 10% of the step, 0 for one that never passes y1
    if max(y) < y0 + 0.9*(y1 - y0)
        s = 'n/a';
    else
        s = sprintf('%.2f', max(0, max(y) - y1)/denom*100);
    end
end


function s = poleStr(p)
    s = '';
    for z = p.'
        if abs(imag(z)) < 1e-9*max(1, abs(real(z)))
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
    late = find(idx & abs(y - target) > tol*target, 1, 'last');
    if isempty(late), rec = 0; else, rec = t(late) - win(1); end
end
