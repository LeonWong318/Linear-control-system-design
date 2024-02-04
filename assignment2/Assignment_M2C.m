clc;
close all;
clear;

% give parameters value
%R = 1;
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
% define B C D matrics 
B = [0,0;0,0;0,1/B_f;K_T/(J_1*R),0;0,0];
% for two cases of different choices of y(t)
C1 = [0,1,0,0,0;0,0,0,0,1];
C2 = [0,0,0,-K_E/R,0;0,D_2/B_f,-D_2/B_f,0,0];
D1 = [0,0;0,0];
D2 = [1/R,0;0,1/B_f];
% question c

W_cc = ctrb (A, B);
W_o1c = obsv(A, C1);
W_o2c = obsv(A, C2);
svd1 = svd(W_cc);
svd2 = svd(W_o1c);
svd3 = svd(W_o2c);
k1 = max(svd1(:))/ (min(svd1(:)));
k2 = max(svd2(:))/ (min(svd2(:)));
k3 = max(svd3(:))/ (min(svd3(:)));
% PBH check the detectivity and controllability
eigenvalue = eig(A);
I = eye(size(A));
PBH1 = [eigenvalue(1)*I - A; C1];
PBH2 = [eigenvalue(1)*I - A; C2];
rref(PBH1);
rref(PBH2);
PBH3 = [eigenvalue(1)*I - A, B];
rref(PBH3);


rank(W_o2c)
rank(W_o1c)

rref(W_o2c)
rref(W_o1c)