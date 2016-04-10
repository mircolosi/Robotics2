function [ c, C ] = compute_cent_coriol_terms( B, q, q_dot )
%COMPUTE_CENT_CORIOL_TERMS Summary of this function goes here
%   Detailed explanation goes here

dim = length(q);
for i = 1:dim
    C(:,:,i) = 1/2*(jacobian(B(:,i),q)+jacobian(B(:,i),q)'-diff(B,q(i)));
    c(i) = q_dot'*C(:,:,i)*q_dot;
end

C = simplify(C);
c = simplify(c);
c = c';

end

