clc; clear all; close all;


%% PARAMETRI

syms q1 q2 q3 q1_dot q2_dot q3_dot pi l l1 l2 l3 dc1 dc2 dc3 m m1 m2 m3 I1 I2 I3 g0 k real
% 
% q = [q1 q2 ]';
% q_dot = [q1_dot q2_dot]';
% 
% pc1=[0;0;0];
% pc2=[(q2-d2)*cos(q1);(q2-d2)*sin(q1);0];
% 
% thetac1=[0;0;q1];
% thetac2=[0;0;q1];
% 
% pc = [pc1 pc2]';
% theta = [thetac1, thetac2]';
% 
% m = [m1 m2]';
% I = [I1 I2]';
% dim = length(q);
% 
% %vettore di gravit√† (es. [0;0;-g0])
% g=[0 -g0 0];



q=[q1 q2 q3]';
q_dot=[q1_dot q2_dot q3_dot]';

%pos dei centri di massa (espressi nelle variabili q1 q2 q3) rispetto a RF0
pc1=[dc1*cos(q1); dc1*sin(q1); 0];
pc2=[l*cos(q1)+dc2*cos(q1+q2);l*sin(q1)+dc2*sin(q1+q2);0];
pc3=[l*cos(q1)+l*cos(q1+q2)+dc3*cos(q1+q2+q3);l*sin(q1)+l*sin(q1+q2)+dc3*sin(q1+q2+q3);0];

% I VETTORI SONO GIA' VERTICALI!!!
pc = [pc1 pc2 pc3];

%angolo attorno a cui ruota il centro di massa del link i. es:
%tethac1=[0;0;q1] significa che c1 ruota di q1 gradi intorno all'asse z del
%RF0
thetac1=[0;0;q1];
thetac2=[0;0;q1+q2];
thetac3=[0;0;q1+q2+q3];

theta = [thetac1 thetac2 thetac3];

% masses and Inertias

m = [m1 m2 m3]';
I = [I1 I2 I3]';

%vettore di gravit‡ (es. [0;0;-g0])
g=[0;-g0;0];


%%

% vc = sym(zeros(3,dim));
% w = sym(zeros(3, dim));
% T = sym(zeros(3,dim));


for i = 1:dim
    vc(:,i) = jacobian(pc(:,i), q)*q_dot;
    w(:,i) = jacobian(theta(:,i), q)*q_dot;
    T(i) = 1/2*m(i)*sqr_norm(vc(:,i))+ 1/2*I(i)*sqr_norm(w(:,i));
end

T = sum(T);

T = collect(simplify(T),q_dot);

B = compute_B(T,q_dot);

[c, C] = compute_cent_coriol_terms(B, q, q_dot);

%% calcola U e g

for i = 1:length(m)
    U(i) = -m(i)*g'*pc(:,i);
end

U=collect(simplify(sum(U)),q)

G=collect(jacobian(U,q)',[g0 k])