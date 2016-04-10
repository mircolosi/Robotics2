clc; clear all; close all;

syms dc1 dc2 dc3 q1 q2 q3 q1_dot q2_dot q3_dot m1 m2 m3 I1 I2 I3 real

vc1 = [-q1_dot*dc1*sin(q1); q1_dot*dc1*cos(q1); 0]; 
w1 = [0; 0; q1_dot];

vc2 = [q2_dot*cos(q1)-q2*q1_dot*sin(q1);q2_dot*sin(q1)+q2*q1_dot*cos(q1);0];
w2 = [0; 0; q1_dot];

vc3 = [ q2_dot*cos(q1)-(q2+dc2)*sin(q1)*q1_dot-dc3*sin(q1+q3)*(q1_dot+q3_dot);
        q2_dot*sin(q1)+(q2+dc2)*cos(q1)*q1_dot+dc3*cos(q1+q3)*(q1_dot+q3_dot);
        0];
w3 = [0; 0; q1_dot+q3_dot];

q = [q1 q2 q3]';
q_dot = [q1_dot q2_dot q3_dot]';
T1 = 1/2*m1*sqr_norm(vc1)+1/2*I1*sqr_norm(w1); 
T2 = 1/2*m2*sqr_norm(vc2)+1/2*I2*sqr_norm(w2);
T3 = 1/2*m3*sqr_norm(vc3)+1/2*I3*sqr_norm(w3);

T_tot = T1+T2+T3;
T_tot = simplify(T_tot)

B = compute_B(T_tot, q_dot);

[c, ~] = compute_cent_coriol_terms(B, q, q_dot);

disp('c ='); disp(c');