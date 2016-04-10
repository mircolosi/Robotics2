function [ v ] = sqr_norm( v, M )
%SQR_NORM Returns the squared norm of a vector v    
    v = v'*M*v;
end

