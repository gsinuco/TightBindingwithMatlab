%Tight binding in a trap and with interactions
function [H_TB] = NLTightBinding_Hamiltonian_driven(t,psi_0)
 
    global N pi J U epsilon delta alpha beta gamma chi omega phi K;

    H_TB = zeros(N,N);

    %Near-Neighbour tunneling matrix
    H_NN = zeros(N,N);
    H_NN = NN_tunneling(K,J);


    % Hubbard mean-field on-site interaction
    H_U = zeros(N,N);
    H_U = onsite(U,psi_0);

    % On-site energy shift due to the trap
    % delt E_i= epsilon i^2 + gamma
    H_U2 = zeros(N,N);
    H_U2 = onsite_shift(epsilon,delta);

    % Site-dependent correction to the tunneling
    % DELTA_I = alpha*i^2 + beta*i + gamma
    H_NN2 = zeros(N,N);
    H_NN2 = NN_tunneling_site(alpha,beta,gamma);

    % density-dependent tunneling rate
    % chi .....
    
    H_NLNNT = zeros(N,N);
    %H_NLNNT = NL_tunneling(chi,psi_0);

    H_DR = zeros(N,N);
    H_DR = NN_Driven_tunneling(alpha,beta,gamma,omega,phi,K,J,t);

    Damping = zeros(N,N);
    Damping = Damping_function(psi_0);
    
    H_TB = zeros(N,N);
    H_TB = H_NN + H_U + H_NN2 + + H_U2 + H_NLNNT + H_DR;% + 1i*0.1*Damping;
    
    
