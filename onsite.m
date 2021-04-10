%% Hubbard onsite interaction

function [H] = onsite(U,psi)

global L N pi;

H = 0.0;

for i=1:N;
    H(i,i) = U*abs(psi(i))^2;
end

