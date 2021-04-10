clc
global L N pi J U epsilon delta alpha beta gamma chi;

pi = 4.0*atan(1.0);
L  = 128;            % number of lattice sites
N  = 10e4;          % Number of particles, for normalization of the wavefunction

J =  -10;
U =   10;

epsilon = 0.0;
delta   = 0.0;

alpha = 0.0;%1.0e-1;
beta  = 0.0;%-alpha;
gamma = 0.0;

chi   = 0.0;

x  = linspace(-L/2,L/2,L);
psi = zeros(L);
%psi = transpose(psi);
%psi = (psi+1);

H = NLTightBinding_Hamiltonian(psi);
Translation=translation(); % Definition of the translation operator: T|i> = |i+1>

%initial state
[V,D]   = eig(H);
psi     = V(:,1);
psi_q   = fftshift(fft(psi));
psi_q   = psi_q/sqrt(sum(abs(psi_q).^2));
psi_q_2 = abs(psi_q).^2;

figure(2)
subplot(1,2,1); 
plot(x,abs(psi).^2,x,abs(V(:,2)).^2,x,abs(V(:,3)).^2);
title('Ground, First, Second: Real space');
xlabel('Position');
ylabel('psi(x)');
drawnow; 

subplot(1,2,2);
plot(x,psi_q_2,'.-');
title('Ground: Momentum space');
xlabel('Momentum');
ylabel('psi(q)');
drawnow;
%%
%%%
%******************************************* 
% ode45 solver settings
%******************************************* 
options=odeset('RelTol',3e-7,'AbsTol',1e-6,'Stats','off','MaxStep',1e-3);


%%
%****************************************************
% Definition of a stable bright soliton (gap soliton)
% Using a selfonsistent solution of the Non-linear
% Schrodinger equation.
%****************************************************
w   = 4;
psi = (1.0/(sqrt(2.0*w)))*sech(x/w);
psi = psi/sqrt(sum(abs(psi).^2));  % normalize
psi = transpose(psi);

N_max    = 32;   % Maximum number of iterations
tol      = 1e-6; % Stop when |psi_n-psi_n-1|<= tol
[psi,mu] = NLSE_SelfconsistenSol(psi,N_max,tol);

figure(4)
plot(x,abs(psi).^2)   
xlabel('Position'); 
ylabel('Density');
%%
%**********************************************
% Construction of the Bogoliubov-de Gennes 
% spectrum
%**********************************************
H_BdG = BdG_Matrix(psi,mu);
[V,D] =eig(H_BdG);
figure(5)
plot(linspace(1,2*L,2*L),real(diag(D)),linspace(1,2*L,2*L),imag(diag(D)))   
%plot(imag(diag(D)))   
xlabel('Index'); 
ylabel('Eigenvalue');



%%

%******************************************* 
% Find the closest stable state
%******************************************* 
Steps  = 128;
tSteps = linspace(0.0,10.0,Steps);

%****************************************************
% Definition of a stable bright soliton (gap soliton)
%****************************************************
w   = 4;
psi = (1.0/(sqrt(2.0*w)))*sech(x/w);
psi = psi/sqrt(sum(abs(psi).^2));  % normalize

psi_x_2 = abs(psi).^2;  % Density psi_x                        
delta_psi = 10;         % Quantifier of the change of the psi
counter = 0;            % iterations counter
dt_gs   = 1e-2;         % time-step to find ground state

while (counter<128)
    counter = counter+1;
    %Imaginary propagation of the Schrodinger equation. 
    [TT,PSI_out] = ode45(@(t,y) DGLFUNC(L,t,y,U,J,epsilon,delta,0),[0 dt_gs/2 dt_gs],psi,options);
    
    psi = PSI_out(3,:);                % take the last time step
    psi = psi/sqrt(sum(abs(psi).^2));  % normalize
    psi_x_3 = abs(psi).^2;             % Density psi_x                       
    if mod(counter,2)==0
        figure(1)
        plot(x,psi_x_3);
        title('Ground state');
        xlabel('Position');
        ylabel('psi(x)');
        drawnow;
    end
    delta_psi = sum(abs(psi_x_3-psi_x_2))              % Has the Wave function changed?
    psi_x_2 = psi_x_3;
end

figure(3)
plot(x,abs(psi).^2,x,abs(V(:,2)).^2,x,abs(V(:,3)).^2)   
xlabel('Position'); 
ylabel('Density');

%%
counter = 0;
dt      = 1.0e-3;
%psi     = V(:,1)+0.1*V(:,2)+0.1*V(:,3);
psi     = psi/sqrt(sum(abs(psi).^2));   % normalize
while (counter < 1280)
    counter = counter+1
    [TT,PSI_out] = ode45(@(t,y) DGLFUNC(L,t,y,U,J,epsilon,delta,1),[0 dt/2 dt],psi,options);
    
    psi = PSI_out(3,:);                    % take the last time step
    %psi = psi/sqrt(sum(abs(psi).^2));  % normalize
    psi_x_3 = abs(psi).^2;                 % Density psi_x                       
    if mod(counter,2)==0
        psi_q = fftshift(fft(psi));
        psi_q = psi_q/sqrt(sum(abs(psi_q).^2));
        psi_q_2 = abs(psi_q).^2;

        figure(3)
        subplot(1,2,1); 
        plot(x,psi_x_3,x,abs(V(:,2)).^2,x,abs(V(:,3)).^2);
        title('Time Evolution');
        xlabel('Position');
        ylabel('psi(x)');
        drawnow; 

        subplot(1,2,2);
        plot(x,psi_q_2,'.-');
        title('Time Evolution');
        xlabel('Momentum');
        ylabel('psi(q)');
        drawnow;
    end
    %delta_psi = sum(abs(psi_x_3-psi_x_2))              % Has the Wave function changed?
    %psi_x_2 = psi_x_3;
end

%%
%[t,psi_t] = ode45(@(t,y) DGLFUNC(L,t,y,U,J,epsilon,delta),tSteps,psi,options);

psi_x = psi_t(end,:);
psi_x = psi_x/sqrt(sum(abs(psi_x).^2));
psi_x_2 = abs(psi_x).^2;                          % Density psi_x    
    
psi_q = fftshift(fft(psi_x));
psi_q = psi_q/sqrt(sum(abs(psi_q).^2));
psi_q_2 = abs(psi_q).^2;
    
%%    
figure(2)
subplot(2,1,1)
plot(x,abs(psi_t).^2)   
xlabel('Position'); 
ylabel('Density');

subplot(2,1,2)
plot(x,psi_q_2)   
xlabel('Momentum'); 
ylabel('Density');


%% 
function dpsi = DGLFUNC(L,t,psi,U,J,epsilon,delta,RealImag)
%but here I need to define H, as I did it above%%
    H_NN = zeros(L,L);
    H_NN = NN_tunneling(J);

    H_U = zeros(L,L);
    H_U = onsite(U,psi);

    H_U2 = zeros(L,L);
    H_U2 = onsite_shift(epsilon,delta);
    
    if(RealImag == 0) % imaginary time evoutlin
        H    =  NLTightBinding_Hamiltonian(psi);
        dpsi = -H*psi;
    end
    if(RealImag == 1) % Real time evolution
        H    = NLTightBinding_Hamiltonian(psi); 
        dpsi = -1i*H*psi;
    end
end

