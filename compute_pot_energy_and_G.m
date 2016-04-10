function [ U, G ] = compute_pot_energy_and_G(pc,q,g,m,k)
%COMPUTE_POT_ENERGY_AND_G Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(m)
    U(i) = -m(i)*g'*pc(:,i);
end
U=collect(simplify(sum(U)),q);
G=collect(jacobian(U,q)',[sym('g0') k]);
end

