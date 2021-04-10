%**************************************************************************
% Here we explore the time evolution of the ground state of the
% tight-binding model (with on-site interactions), after a sudden change 
% of the sign of the tunneling rate.
%**************************************************************************

clc
global L N pi J U epsilon delta alpha beta gamma chi omega phi K;

pi= 4.0*atan(1.0);
%L = 17;            % number of lattice sites
N = 256;%10e4;          % Number of particles, for normalization of the wavefunction

%Hubbard model parameters
J = 0.0227;%10.0;
U =  0.0;%0.5;

%Harmonic trap: local energy shifts
% deltaE_i = epsilon*x_i^2+ delta
epsilon = 2.5e-5;%0.125;
delta   = 0.0402;%0.0;

%Harmonic trap: local tunneling rate
% DELTA_i = alpha*x_i^2 + beta*x_i + gamma
alpha = 0.0;%1.0e-1;
beta  = 0.0;%-alpha;
gamma = 0.0;

%density dependent tunnneling:

chi   = 0.0;

%Shaking parameters
% cos(k*(x+A*cos(omega*t + phi)))
omega = 1e3;
phi   = 0.0;

% Renormalization argument
% K = A/m omega
K     = 0.0;%1.74;


x   = linspace(-N/2,N/2,N);
psi = zeros(N);

H = NLTightBinding_Hamiltonian(psi);
Translation=translation(); % Definition of the translation operator: T|i> = |i+1>

%initial state
[V,D]   = eig(H);
psi     = V(:,1);
psi_q   = fftshift(fft(psi));
psi_q   = psi_q/sqrt(sum(abs(psi_q).^2));
psi_q_2 = abs(psi_q).^2;

figure(1)
    subplot(1,2,1); 
    plot(x,abs(psi).^2,x,abs(V(:,2)).^2,x,abs(V(:,3)).^2);
    axis([-N/2 N/2 0 inf]); 
    title('Ground, First, Second: Real space');
    xlabel('Position');
    ylabel('psi(x)');
    drawnow; 

    subplot(1,2,2);
    plot(x,psi_q_2,'.-');
    axis([-N/2 N/2 0 inf]);     
    title('Ground: Momentum space');
    xlabel('Momentum');
    ylabel('psi(q)');
    drawnow;

%%%******************************************* 
% ode45 solver settings
%******************************************* 
options=odeset('RelTol',3e-7,'AbsTol',1e-6,'Stats','off','MaxStep',1e-3);

%****************************************************
% Definition of a stable ground state
% Using an imaginary-time propagation of the Non-linear
% Schrodinger equation.
%****************************************************

Steps  = 128;
tSteps = linspace(0.0,10.0,Steps);

psi_x_2   = abs(psi).^2;  % Density psi_x                        
delta_psi = 10;         % Quantifier of the change of the psi
counter   = 0;            % iterations counter
dt_gs     = 1.0e-1;         % time-step to find ground state

while (delta_psi > 1e-6)
    counter = counter+1
    %Imaginary propagation of the Schrodinger equation. 
    [TT,PSI_out] = ode45(@(t,y) DGLFUNC(N,t,y,0),[0 dt_gs/2 dt_gs],psi,options);
    
    psi = PSI_out(3,:);                       % take the last time step
    psi = psi/sqrt(sum(abs(psi).^2));         % normalize
    mu  = -log(1.0/sqrt(sum(abs(PSI_out(3,:)).^2)))/dt_gs; % Chemical potential
    psi_x_3 = abs(psi).^2;                    % Density psi_x                       
    if mod(counter,10)==0
        figure(2)
        plot(x,psi_x_3);
        hold on
        title('Ground state');
        xlabel('Position');
        ylabel('psi(x)');
        drawnow;
    end
    delta_psi = sum(abs(psi_x_3-psi_x_2)); % Has the Wave function changed?
    psi_x_2 = psi_x_3;
end

psi_0 = psi
psi_q   = fftshift(fft(psi));
psi_q   = psi_q/sqrt(sum(abs(psi_q).^2));
psi_q_2 = abs(psi_q).^2;

figure(3)
    subplot(1,2,1); 
    plot(x,abs(psi_0).^2,x,abs(V(:,2)).^2,x,abs(V(:,3)).^2)   
    axis([-N/2 N/2 0 inf]); 
    title('Ground, First, Second: Real space');
    xlabel('Position');
    ylabel('psi(x)');
    drawnow; 

    subplot(1,2,2);
    plot(x,psi_q_2,'.-');
    axis([-N/2 N/2 0 inf]); 
    title('Ground: Momentum space');
    xlabel('Momentum');
    ylabel('psi(q)');
    drawnow;

%**********************************************
% Construction of the Bogoliubov-de Gennes 
% spectrum
%**********************************************
H_BdG = BdG_Matrix(psi,mu);
[V,D] =eig(H_BdG);
figure(4)
plot(linspace(1,2*N,2*N),real(sort(diag(D))),'.-',linspace(1,2*N,2*N),imag(diag(D)),'.-')   
%plot(imag(diag(D)))   
xlabel('Index'); 
ylabel('Eigenvalue');


%%
%%%******************************************* 
% ode45 solver settings
%******************************************* 
options=odeset('RelTol',1e-7,'AbsTol',1e-7,'Stats','off','MaxStep',1.0e-2);

K = 3.26;%4.2;
%J = 10;

counter = 0;
dt      = 1.0e-3;

%phase = exp(0.001*1i*2*pi*rand(1,L));
%psi   = psi.*phase;
%%
psi = psi_0;
psi_q = fftshift(fft(psi));
psi_q = psi_q/sqrt(sum(abs(psi_q).^2));
psi_q_2 = abs(psi_q).^2;

figure(8)
    subplot(1,2,1); 
    plot(x,abs(psi),'.-')%,x,abs(V(:,2)).^2,x,abs(V(:,3)).^2);
    axis([-N/2 N/2 0 0.55]); 
    title('Time Evolution');
    xlabel('Position');
    ylabel('psi(x)');
    drawnow; 

    subplot(1,2,2);
    plot(x,psi_q_2,'.-');
    axis([-N/2 N/2 0 0.55]); 
    title('Time Evolution');
    xlabel('Momentum');
    ylabel('psi(q)');
    drawnow;

%%
counter = 0

Steps  = 256;
tSteps = linspace(0,1,Steps);
% buffer
psi_x_vect   = zeros( 2*length(tSteps)-1, N );
psi_x_2_vect = psi_x_vect;
psi_q_vect   = psi_x_vect;
psi_q_2_vect = psi_x_vect;
%psi_exp_2_vect = zeros( 2*length(tSteps)-1, length(psi_exp_2) );
 
while (counter < 10)
    counter = counter+1
    
    %[TT,PSI_out] = ode45(@(t,y) DGLFUNC(N,t,y,1),[0 1],psi,options);
    [TT,PSI_out] = ode45(@(t,y) DGLFUNC(N,t,y,1),tSteps,psi,options);
    
    for j1=1:length(tSteps)
        psi_x_vect(j1,:) =  PSI_out(j1,:);
        psi_x_vect(j1,:) =  psi_x_vect(j1,:)/sqrt(sum(abs(psi_x_vect(j1,:)).^2));
        psi_x_2_vect(j1,:) = abs( psi_x_vect(j1,:) ).^2;  
        
        psi_q_vect(j1,:) = fftshift(fft(psi_x_vect(j1,:)));
        psi_q_vect(j1,:) = psi_q_vect(j1,:)/sqrt(sum(abs(psi_q_vect(j1,:)).^2));
        psi_q_2_vect(j1,:) = abs(psi_q_vect(j1,:)).^2;        
            figure(7)
            subplot(1,2,1); 
            plot(x,psi_x_2_vect(j1,:))%,x,abs(V(:,2)).^2,x,abs(V(:,3)).^2);
            axis([-N/2 N/2 0 0.55]); 
            title('Time Evolution');
            xlabel('Position');
            ylabel('psi(x)');
            drawnow; 

            subplot(1,2,2);
            plot(x,psi_q_2_vect(j1,:),'.-');
            axis([-N/2 N/2 0 0.55]); 
            title('Time Evolution');
            xlabel('Momentum');
            ylabel('psi(q)');
            drawnow;
        %psi_exp_2_vect(j1,:) = conv(PSI_Exp_Start,psi_q_2_vect(j1,:))';
        %psi_exp_2_vect(j1,:) = psi_exp_2_vect(j1,:)/sum(psi_exp_2_vect(j1,:));
    end    
    
    psi = PSI_out(end,:);                 % take the last time step
    %psi = psi/sqrt(sum(abs(psi).^2));  % normalize
    psi_x_3 = abs(psi).^2;              % Density psi_x                       
    if mod(counter,1)==0
        psi_q = fftshift(fft(psi));
        psi_q = psi_q/sqrt(sum(abs(psi_q).^2));
        psi_q_2 = abs(psi_q).^2;

        figure(7)
            subplot(1,2,1); 
            plot(x,psi_x_3)%,x,abs(V(:,2)).^2,x,abs(V(:,3)).^2);
            axis([-N/2 N/2 0 0.55]); 
            title('Time Evolution');
            xlabel('Position');
            ylabel('psi(x)');
            drawnow; 

            subplot(1,2,2);
            plot(x,psi_q_2,'.-');
            axis([-N/2 N/2 0 0.55]); 
            title('Time Evolution');
            xlabel('Momentum');
            ylabel('psi(q)');
            drawnow;
    end
    %delta_psi = sum(abs(psi_x_3-psi_x_2))              % Has the Wave function changed?
    %psi_x_2 = psi_x_3;
end


%% 
function dpsi = DGLFUNC(N,t,psi,RealImag)
     
    if(RealImag == 0) % imaginary time evoutlin
        H    =  NLTightBinding_Hamiltonian(psi);
        dpsi =  -H*psi;
    end
    if(RealImag == 1) % Real time evolution
        %H    =     NLTightBinding_Hamiltonian(psi); 
        H    =     NLTightBinding_Hamiltonian_driven(t,psi); 
        dpsi = -1i*H*psi;
    end
    
end

