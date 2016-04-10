function [ T_tot, w, v, vc ] = compute_kin_energy( R, r, rc, m, I, q_dot, type, n )
%COMPUTE_POTENTIAL_ENERGY Summary of this function goes here
%   Detailed explanation goes here

w(:,1) = sym('[0,0,0]')';
v(:,1) = sym('[0,0,0]')';
vc(:,1) = sym('[0,0,0]')';

for i=2:n+1
    
    if strcmp(type(i-1), 'r')
        sigma = 0;
    else
        sigma = 1;
    end
    w_iminus1_i = (w(:,i-1)+(1-sigma)*q_dot(i-1)*[0 0 1]');
    w(:,i) = R(:,:,i-1)'* w_iminus1_i;
    w(:,i) = simplify(w(:,i));
    v(:,i) = R(:,:,i-1)'* (v(:,i-1)+sigma*q_dot(i-1)*[0 0 1]'+cross(w_iminus1_i, r(:,i-1)));
    v(:,i) = simplify(v(:,i));
    vc(:,i) = v(:,i) + cross(w(:,i), rc(:,i-1));
    vc(:,i) = simplify(vc(:,i));
    T(i-1) = 1/2*m(i-1)*sqr_norm(vc(:,i))+1/2*I(i-1)*sqr_norm(w(:,i));
    T(i-1) = simplify(T(i-1));
end

disp('w = ');disp(w);
disp('v = ');disp(v);
disp('vc = ');disp(vc);

disp('T = ');pretty(T);

T_tot = simplify(sum(T));

end

