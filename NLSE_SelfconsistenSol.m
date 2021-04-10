%NL Schordinger self-consistent solution
function [psi_sc,mu] = NLSE_SelfconsistenSol(psi_0,N_max,tol)


    global L N pi J U epsilon delta alpha beta gamma chi;

    counter = 0;
    psi_sc = zeros(L);
    x  = linspace(-64,64,L);

    while (counter < N_max)% | delta_psi > tol)

        H         = NLTightBinding_Hamiltonian(psi_0);
        [V,D]     = eig(H);
        [M,I]     = max(abs(transpose(V)*psi_0));
        psi_sc    = V(:,I);
        delta_psi = sum(abs(abs(psi_sc).^2-abs(psi_0).^2)); % Has the Wave function changed?
        counter   = counter + 1
        psi_0     = psi_sc;
        mu        = D(I,I);
        figure(3)
        plot(x,abs(psi_sc).^2)   
        hold on
        xlabel('Position'); 
        ylabel('Density');
    end
    
    
     
    