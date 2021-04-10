%% Hubbard onsite interaction shift due to the presence of the trap

function [H] = onsite_shift(epsilon,delta)

global L N pi;

H = 0.0;

for i=1:N;
    site_i = -(N-1)/2 + (i-1);
    H(i,i) = epsilon*site_i*site_i + delta;
end
