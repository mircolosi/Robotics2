function [ B ] = compute_B( T_tot, q_dot )
%COMPUTE_B Summary of this function goes here
%   Detailed explanation goes here

[coeff, terms] = coeffs(T_tot, q_dot);

coeff = simplify(coeff');
coeff = collect(coeff, 1/2);

terms = terms';
clen = length(coeff);
dim = length(q_dot);
B=sym(zeros(max(size(q_dot))));
for l = 1:clen
    
    lterm=char(terms(l));
    i_indexes=strfind(lterm,'q');
    if size(i_indexes, 2)>1
        %index 1
        i=str2num(lterm(i_indexes(1)+1));
        %index 2
        j=str2num(lterm(i_indexes(2)+1));
        B(i,j) = coeff(l);
        B(j,i) = coeff(l);
    else
        %index 1
        i=str2num(lterm(i_indexes(1)+1));
        B(i,i) = 2*coeff(l);   
    end
end

B = simplify(B);

end

