%% Assignment M03
clc;
close all;
clear;

K_E = 10^-1;
K_T = 10^-1;
J_1 = 10^-5;
J_2 = 4 * 10^-5;
B_f = 2 * 10^-3;
D_2 = 2;
D_1 = 20;
R = 1;

% define A matrics

A = [0,0,0,1,0;
    0,0,0,0,1;
    0,D_2/B_f,-D_2/B_f,0,0;
    -D_1/J_1,D_1/J_1,0,-(K_E*K_T)/(J_1*R),0;
    D_1/J_2,-(D_1+D_2)/J_2,D_2/J_2,0,0];

B = [0,0;0,0;0,1/B_f;K_T/(J_1*R),0;0,0];
C1 = [0,1,0,0,0;0,0,0,0,1];
C2 = [0,0,0,-K_E/R,0;0,D_2/B_f,-D_2/B_f,0,0];
% calculate matrix Ad and Bd
Ts = 0.001;
Ad = expm(A*Ts);
syms s t;
variable = inv (s*eye(size(A))-A);
exp_At = vpa(ilaplace(variable));
Bd = int(exp_At, t, 0, 1e-3) * B;
doubleBd = double(Bd)

%% Question a)

N = doubleBd;

% Calculate noise variaces
sigma2_va   = calcVariance(0.3, 0.997);
sigma2_Te   = calcVariance(0.1, 0.997);
sigma2_phi2 = calcVariance(0.02, 0.997);
sigma2_w2   = calcVariance(0.01, 0.997);
% Question a/b)
R = diag([sigma2_va sigma2_Te sigma2_phi2 sigma2_w2])


%% Question c)

% Calculate Kalman gain
sysmodel = ss(Ad, [doubleBd N], C1, 0, Ts)
Qm = diag([sigma2_va sigma2_Te]); % input noise
Rm = diag([sigma2_phi2 sigma2_w2]); % output noise
% kest, kalman gain and covariance value
[kest,L_kalman,P] = kalman(sysmodel, Qm , Rm)

% Check observer stability (observer eigenvalues)
eig(Ad-L_kalman*C1)

%% Question d)

x0 = [0 0 0 0 0]';

% C1_I = [1 1];
% C1_fb = C1_I*C1;
C1_fb = [0 0 0 0 1];

% Set extended system model
Ad_edirect = [Ad  0*Ad*C1_fb';
       -C1_fb eye(size(C1_fb,1))]
Bd_edirect = double([Bd;
        0*(C1_fb*Bd)])
% C1_e= [C1_fb 0*(C1_fb*C1_fb')]
C1_e= [0 0 0 0 1 0];
Ae = [A  0*A*C1_fb';
    -C1_e];
Be = [B;
    0*(C1_fb*B)];

% Bd_edouble = double(Bd);
% sysmodel_e = ss(Ad_e,Bd_edouble,C1_e, 0, Ts);
sys = ss(Ae, Be, C1_e, 0);
sysdt = c2d(sys, Ts);
Ad_e = double(sysdt.A);
Bd_e = double(sysdt.B);
svd_pos = svd(ctrb(Ad_e,Bd_e));
cond_number = max(svd_pos)/min(svd_pos)
% Calculate LQI controller
% Qx = diag([0 0 0 10 1000 10]); % set suitable parameters
Qu = diag([100 1000]);
% Qx = diag([1 1 10 0.001 1 0.001]);
% Qu = diag([2 2])*1e6;
% Be_double = double(Be);
Qx = diag([0 0 0 10 1000 10]); % set suitable parameters
% Qu = diag([1000 0.1]);

[Klqr,S,e]=dlqr(Ad_e, Bd_e,  Qx, Qu)
% [Klqr,S,e] = dlqr(Ad_e, Bd_e, Qx, Qu);
% [Klqr,S,e] = lqrd(Ad_e, Bd_edouble, Qx, Qu, Ts);
% Set variables for Simulink
K_fb  = Klqr(:,1:length(Ad));
K_I = Klqr(:,length(Ad)+1:end)*20;


% function
function sigma2_val = calcVariance(ub, realization_ub)
    syms sigma2 x
    eq = realization_ub == int(1/sqrt(2*pi*sigma2)*exp(-(x^2)/(2*sigma2)),x,-ub,ub);
    sigma2_val = double(solve(eq,sigma2));
    mu = 0;
    pd = makedist('Normal',mu,abs( sqrt(sigma2_val)));
    v_aw = linspace(-ub*2,ub*2,100);
    y = pdf(pd,v_aw);
%     figure('Color','white');
%     plot(v_aw,y);
end