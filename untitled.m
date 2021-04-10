clc
global L N pi J U epsilon delta alpha beta gamma chi omega phi K;

pi= 4.0*atan(1.0);
L = 65;            % number of lattice sites
N = 10e4;          % Number of particles, for normalization of the wavefunction

%Hubbard model parameters
J = 1.0;
U =  0.5;

%Harmonic trap: local energy shifts
% deltaE_i = epsilon*x_i^2+ delta
epsilon = 1.0;
delta   = 0.0;


%%%******************************************* 
% ode45 solver settings
%******************************************* 
options=odeset('RelTol',3e-7,'AbsTol',1e-6,'Stats','off','MaxStep',1e-3);


while (delta_psi > 1e-6)
    counter = counter+1
    %Imaginary propagation of the Schrodinger equation. 
    [TT,PSI_out] = ode45(@(t,y) DGLFUNC(L,t,y,0),[0 dt_gs/2 dt_gs],psi,options);
    
    psi = PSI_out(3,:);                       % take the last time step
    figure(2)
        plot(t,psi(1),t,psi(2));
        hold on
        title('SC');
        xlabel('time');
        ylabel('k,x)');
        drawnow;
    end
end


%% 
function dpsi = DGLFUNC(L,t,psi,RealImag)
   
        dpsi(1) = ;    
        dpsi(2) = ;
end


