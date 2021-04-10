%% Bogoliuvob-deGennes Matrix

function [H_BdG] = BdG_Matrix(psi_0,mu)

    global L N pi J U epsilon delta alpha beta gamma chi;

    H_BdG = zeros(2*N,2*N);
    % Hubbard mean-field on-site interaction
    H_U = zeros(N,N);
    H_U = onsite(U,psi_0);
    
    H = NLTightBinding_Hamiltonian(psi_0);
    H_BdG(  1:N  ,  1:N)   =   H + H_U  - mu*eye(N);
    H_BdG(N+1:2*N,N+1:2*N) = -(H + H_U) + mu*eye(N);
    H_BdG(  1:N,  N+1:2*N) = -U*abs(diag(psi_0)).^2;
    H_BdG(N+1:2*N,  1:N)   =  U*abs(diag(psi_0)).^2;



