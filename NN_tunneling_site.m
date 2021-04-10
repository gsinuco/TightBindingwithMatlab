%% Site-dependent tunnelling interaction

function [H] = NN_tunneling_site(alpha,beta,gamma)

global L N pi;

H = zeros(N,N);


for i=1:N;
    site_i = -(N-1)/2 + (i-1);
    if(i<N)
        H(i,i+1) = alpha*site_i^2 + beta*site_i ;%+ gamma;
    end
    if(i>1)
        H(i,i-1) = alpha*site_i^2 - beta*site_i ;%+ gamma;    
    end
    if(i==1)
        H(i,N) = alpha*site_i^2 - beta*site_i ;%+ gamma;            
        H(i,2) = alpha*site_i^2 + beta*site_i ;%+ gamma;            
    end
    if(i==N)
        H(i,N-1) = alpha*site_i^2 - beta*site_i ;%+ gamma;            
        H(i,1)   = alpha*site_i^2 + beta*site_i ;%+ gamma;            
    end
    %H(i+1,i) = H(i,i+1);
end

