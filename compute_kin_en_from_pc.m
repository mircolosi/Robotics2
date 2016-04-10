function [ T ] = compute_kin_en_from_pc( pc, theta, q, q_dot, m, I )
%COMPUTE_KIN_EN_FROM_PC Summary of this function goes here
%   Detailed explanation goes here

dim = length(q);

for i = 1:dim
    vc(:,i) = jacobian(pc(:,i), q)*q_dot;
    w(:,i) = jacobian(theta(:,i), q)*q_dot;
    T(i) = 1/2*m(i)*sqr_norm(vc(:,i))+ 1/2*I(i)*sqr_norm(w(:,i));
end

T = sum(T);

T = collect(simplify(T),q_dot);

end

