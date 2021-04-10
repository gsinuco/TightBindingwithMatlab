%% Site and density dependent tunneling interaction

function [H] = NL_tunneling_site(chi,psi)

global L N pi;

H = 0.0;

for i=2:N-1;
    site_i = -(N-1)/2 + (i-1);
    H(i,i+1) = chi*psi(i)*psi(i+1);
    H(i+1,i) = chi*psi(i)*psi(i+1);
end
H(1,2)   = chi;
H(2,1)   = chi;
%H(L,L-1) = chi;
%H(L-1,L) = chi;