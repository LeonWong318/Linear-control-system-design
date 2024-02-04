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
% question d
Ts = 0.001;
Ad = expm(A*Ts)

% question e
syms s t;
variable = inv (s*eye(size(A))-A);
exp_At = vpa(ilaplace(variable));
Bd = int(exp_At, t, 0, 1e-3) * B
double(Bd)

% question f

%Singular value decomposition of controllability matrix (discrete time)
svd1 = svd(ctrb(Ad, Bd))
%Singular value decomposition of observability matrix case 1(discrete time)
svd2 = svd(obsv(Ad, C1))
%Singular value decomposition of observability matrix case 2(discrete time)
svd3 = svd(obsv(Ad, C2))
% stability check
eigenvalue2 = abs(eig(Ad))