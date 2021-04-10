%% Near_neighbour tunneling

function [H] = NN_Driven_tunneling(alpha,beta,gamma,omega,phi,K,J,t)

    global L N;

    H = zeros(N,N);

    for i=1:N;
        site_i = -(N-1)/2 + (i-1);
        if(i<N)
            J_i = J + alpha*site_i^2 + beta*site_i + gamma;        
            aux = 0.0;
            for n=-3:3            
                if n ~= 0
                    aux = aux - J_i*exp(-1i*n*(omega*t+phi))*besselj(n,K);
                end
            end
            H(i,i+1) = aux;
        end
        if(i>1)
            J_i = J + alpha*site_i^2 - beta*site_i + gamma;        
            aux = 0.0;
            for n=-3:3            
                if n ~= 0
                    aux = aux - J_i*exp(1i*n*(omega*t+phi))*besselj(n,K);
                end
            end
            H(i,i-1) = aux;
        end
        if(i==1)
            J_i    = alpha*site_i^2 - beta*site_i + gamma;            
            aux = 0.0;
            for n=-3:3            
                if n ~= 0
                    aux = aux - J_i*exp(1i*n*(omega*t+phi))*besselj(n,K);
                end
            end
            H(1,N) = aux;
        end
        if(i==N)
            J_i       = alpha*site_i^2 + beta*site_i + gamma;            
            aux = 0.0;
            for n=-3:3            
                if n ~= 0
                    aux = aux - J_i*exp(-1i*n*(omega*t+phi))*besselj(n,K);
                end
            end
            H(N,1) = aux;
        end
    end



