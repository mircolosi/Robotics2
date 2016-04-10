function [ A, T0E ] = DH_Matrix(d, theta, r, alpha)
%DH_MATRIX function returns the A matrix and the resulting matrix of the 
%homogeneus transformation given a of the direct kinematics given the
%parameters.
%   n is the number of joints
%   d is the vector of the distances d along z of the 2 origins of the two
%   consecutive frames
%   theta is the vector of the angles around z between x(i-1) and x(i)
%   r is the radius of the rotation from O(i) to O(i+1)
%   alpha is the twist angle between z(i-1) and z(i) around x(i)
%-------------RACCOMANDAZIONE--------------
%   UTILIZZA sym(pi) al posto di pi (oppure all'inizio dichiara pi=sym(pi))
%   per evitare che il pi venga calcolato e dia problemi di approssimazione
%   UTILIZZA la trasposizione puntata o dichiara i simboli reali 

n = max(size(d));
    for i=1:n
        A(:,:,i)=[  cos(theta(i))   -cos(alpha(i))*sin(theta(i))    sin(alpha(i))*sin(theta(i))     r(i)*cos(theta(i))  ;
                    sin(theta(i))   cos(alpha(i))*cos(theta(i))     -sin(alpha(i))*cos(theta(i))    r(i)*sin(theta(i))  ;
                    0               sin(alpha(i))                   cos(alpha(i))                   d(i)                ;
                    0               0                               0                               1                   ];
    end
    
    T0E=eye(4);
    
    for i=1:n
        T0E=T0E*A(:,:,i);
    end
    
end

