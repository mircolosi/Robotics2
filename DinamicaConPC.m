clc; clear all; close all;

syms q1 q2 q1_dot q2_dot g0 m1 m2 I1 I2 k real

pc1 = [0, 0, 0]';
pc2 = [q2*cos(q1), q2*sin(q1), 0]';
pc = [pc1 pc2];

theta1 = [0 0 q1]';
theta2 = [0 0 q1]';
theta = [theta1 theta2];

g = [0 -g0 0]';

q = [q1 q2]';
q_dot = [q1_dot q2_dot]';

m = [m1 m2]';
I = [I1 I2]';


%%
T = compute_kin_en_from_pc(pc, theta, q, q_dot, m, I);
disp('T= '); disp(T);
B = compute_B(T,q_dot);
disp('B= '); disp(B);
[c, C] = compute_cent_coriol_terms(B, q, q_dot);
disp('c= '); disp(c);

%% calcola U e g

[U, G] = compute_pot_energy_and_G(pc,q,g,m,k);

disp('U ='); disp(U);
disp('G ='); disp(G);