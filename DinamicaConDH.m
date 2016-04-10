clc; clear all; close all;

syms q1 q2 q3 q1_dot q2_dot q3_dot pi l1 l2 l3 d1 d2 d3 m1 m2 m3 I1 I2 I3 real

d           = [0        q2];
theta       = [q1       0];
r           = [0        0];       
a           = [pi/2     0];
type        = ['r'      'p' ];
q           = [q1       q2 ]';
q_dot       = [q1_dot   q2_dot]';

rc          = [ 0   0   0; 
                0   0   -d2]';

m           = [m1 m2 ];
I           = [I1 I2 ];

[A,~] = DH_Matrix(d, theta, r, a);
R = A(1:3,1:3,:);
r = A(1:3,4,:);

n = length(theta);

[T_tot, w, v, vc] = compute_kin_energy(R, r, rc, m, I, q_dot, type, n);
disp('T= '); disp(T_tot);

B = compute_B(T_tot, q_dot);
disp('B= '); disp(B);

[c, ~] = compute_cent_coriol_terms(B, q, q_dot);

disp('c ='); disp(c);