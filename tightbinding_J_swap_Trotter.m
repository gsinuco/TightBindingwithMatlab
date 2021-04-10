%**************************************************************************
% Here we explore the time evolution of the ground state of the
% tight-binding model (with on-site interactions), after a sudden change 
% of the sign of the tunneling rate.
%
% Here we ignore interactions and evaluate the time-evolution operator
% during a period. 
%
% U(t_0+T,t_0) = prod exp(-i H(t) DELTA_t/hbar)
%
%**************************************************************************

clc
global N pi J U epsilon delta alpha beta gamma chi omega phi K;

%******************************************* 
% constants
%******************************************* 
c      = 299792458;           % Speed of light [m/s]
kB     = 1.3806503e-23;       % Boltzmann constant [J/K]
hbar   = 1.05457148e-34;      % Plank constant [Js]
muB    = 9.27400949e-24;      % Bohr magneton [J/T]
g      = 9.80553;             % Acceleration of gravity [m/s^2]
u      = 1.6605e-27;          % Atomic mass unit [kg]
a0     = 0.5291772083e-10;    % Bohr's radius in [m]
pi     = 4.0*atan(1.0);

mCs    = 133*u;               % Mass of 133Cs [kg]
lLat   = 1064.49463e-9;       % Lattice wavelength [m]
kLat   = 2*pi/lLat;           % Lattice wavevektor [1/m]
omgRec = hbar*kLat^2/2/mCs;   % Recoil frequency [Hz]
G      = 2*kLat;              % Brilloin zone width [1/m]
Er     = hbar^2*kLat^2/2/mCs; % Recoil energy [J]
dLat   = lLat/2;              % Lattice site distance [m]


%******************************************* 
% Modulation and time steps
%******************************************* 
modPeriod =  1.0e-3;        % modulation period in s
modOmega  =  2*pi/modPeriod; % modulation angular frequency
nPeriods  = 100;             % Number of periods for simulation
Steps     = 11;              % Saved steps per BO
dt        = modPeriod;       % Single time step dt
nSaved = 1 + nPeriods*Steps; % number saved wave functions

%******************************************* 
% Physical system, general parameters 
%******************************************* 
N       = 65;       % Number of Wells, must be odd
s       = 10;         % Lattice Depth in Er
a       = 0.0*a0;      % Scattering length
omgLong = 2*pi * 5;  % Longitudinal trap frequency
omgPerp = 2*pi * 60; % radial, perpendicular trap freq.

nAtoms  = 1E4;       % Number of atoms
expTime = 40e-3;     % End expansion time

% length scales
sigLong = sqrt(hbar/mCs/omgLong); % Longitudinal oscillator length
sigPerp = sqrt(hbar/mCs/omgPerp); % Radial oscillator length

% assume that vertical beam is strong, large radial frequency
omgMean = (omgPerp*omgPerp*omgLong)^(1/3); % Geometric averaged trap frequency
sigMean = sqrt(hbar/mCs/omgMean);          % Geometric averaged oscillator length
mu      = hbar*omgMean/2*((15*nAtoms*a/sigMean)^(2/5)); % Chemical Potential for TF-approx.
rTFLong = sqrt(2*mu/mCs/omgLong/omgLong);               % Longitudinal TF-radius

omgLatt = 2*omgRec*sqrt(s);       % Harmonic frequency in one well
sigLatt = sqrt(hbar/mCs/omgLatt); % Oscillator length in one well

MaxExpqr = hbar*kLat/mCs*expTime; % Expansion faktor

%%%******************************************* 
% ode45 solver settings
%******************************************* 
options  = odeset('RelTol',3e-7,'AbsTol',1e-6,'Stats','off','MaxStep',1e-3);

%%%******************************************* 
% tight-binding parameters from the optical lattice 
%******************************************* 

% Renormalization argument
K     = 3.26;  % Argument of the Bessel function 
omega = modOmega*modPeriod; % scaled driving frequency
                                  % the unit of time is the driving period
phi   = 0; % phase of the driving      

[J,U,alpha,beta,gamma,epsilon,delta,chi] = OL2TB_parameters(nAtoms , s , lLat , omgPerp, omgLong , mCs , a );


% The Hamiltonian is scaled with Er
% The unit of time is the driving period
% U(t_0+T,t_0) = prod exp(-i H(t) DELTA_t/hbar)
steps  = 256;
U_T    = eye(N);
t      = 0.0;
dt     = 1.0/steps;
psi    = zeros(N);

for i=1,steps;    
    H    = NLTightBinding_Hamiltonian_driven(t,psi); 
    H_dt = H*dt*Er*modPeriod*hbar;
    U_T  = U_T*expm(-1i*H_dt);
    t    = t + dt;
end


[V,D]   = eig(U_T);

%%

% Renormalization argument
K     = 0.0;%3.26;  % Argument of the Bessel function 
omega = modOmega*modPeriod; % scaled driving frequency
                                  % the unit of time is the driving period
phi   = 0; % phase of the driving      

Schrodinger_Scale = modPeriod*Er/hbar; % Scale of the RHS of the Schrodinger equation


x   = linspace(-N/2,N/2,N);
psi = zeros(N);
N

H = NLTightBinding_Hamiltonian(psi);
Translation=translation(); % Definition of the translation operator: T|i> = |i+1>


%%%******************************************* 
%%
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
%%
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
    [TT,PSI_out] = ode45(@(t,y) DGLFUNC(Schrodinger_Scale,N,t,y,0),[0 dt_gs/2 dt_gs],psi,options);
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

psi_0 = psi;
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

K = 3.26;
J = 10;

counter = 0;
dt      = 1.0e-2;

%phase = exp(0.001*1i*2*pi*rand(1,L));
%psi   = psi.*phase;

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

Steps  = 64*16;
tSteps = linspace(0,1,Steps);
% buffer
psi_x_vect   = zeros( 2*length(tSteps)-1, N );
psi_x_2_vect = psi_x_vect;
psi_q_vect   = psi_x_vect;
psi_q_2_vect = psi_x_vect;
%psi_exp_2_vect = zeros( 2*length(tSteps)-1, length(psi_exp_2) );
 
while (counter < 1)
    counter = counter+1
    
    %[TT,PSI_out] = ode45(@(t,y) DGLFUNC(N,t,y,1),[0 1],psi,options);
    [TT,PSI_out] = ode45(@(t,y) DGLFUNC(Schrodinger_Scale,N,t,y,1),tSteps,psi,options);
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
function dpsi = DGLFUNC(Schrodinger_Scale,N,t,psi,RealImag)
     
    if(RealImag == 0) % imaginary time evoutlin
        H    =  NLTightBinding_Hamiltonian(psi);
        dpsi =  -Schrodinger_Scale*H*psi;
    end
    if(RealImag == 1) % Real time evolution
        %H    =     NLTightBinding_Hamiltonian(psi); 
        H    =     NLTightBinding_Hamiltonian_driven(t,psi); 
        dpsi = -1i*Schrodinger_Scale*H*psi;
    end
    
end

