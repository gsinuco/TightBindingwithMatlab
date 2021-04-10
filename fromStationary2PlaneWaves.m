function [U_2PW]=fromStationary2PlaneWaves(i_start)
    global L;
    
    U_2PW = zeros(L,L);
    
    U_2PW(1,1) = 1;
    U_2PW(L,L) = 1;
    if(mod(L,2)==1) 
        i_start = 0;
    end
    if(mod(L,2)==0) 
        i_start = 1;
    end
    for i=i_start:L/2-1
        if(i_start==1)
            U_2PW(2*i,  2*i)   =  1.0/sqrt(2);
            U_2PW(2*i+1,2*i)   =   1i/sqrt(2);
            U_2PW(2*i,  2*i+1) =  1.0/sqrt(2);
            U_2PW(2*i+1,2*i+1) =  -1i/sqrt(2);        
        end
        if(i_start==0)
            U_2PW(2*i+1,  2*i+1)   =  1.0/sqrt(2);
            U_2PW(2*i+2,  2*i+1)   =   1i/sqrt(2);
            U_2PW(2*i+1,  2*i+2)   =  1.0/sqrt(2);
            U_2PW(2*i+2,  2*i+2)   =  -1i/sqrt(2);        
        end
         
%        V_new(:,2*i)   = V(:,2*i)+1i*V(:,2*i+1);
%        V_new(:,2*i+1) = V(:,2*i)-1i*V(:,2*i+1);
    end 

    
    

end