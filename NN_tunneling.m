%% Near_neighbour tunneling

function [H] = NN_tunneling(K,J)

global L N pi;

H = 0.0;

for i=1:N-1;
    H(i,i+1) = -J*besselj(0,K);
    H(i+1,i) = -J*besselj(0,K);
end

% Defining periodic Boundary conditions
H(N,1) = -J*besselj(0,K);
H(1,N) = -J*besselj(0,K);

