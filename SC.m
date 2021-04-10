clc
global L N pi J U epsilon delta alpha beta gamma chi omega phi K;

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
modPeriod =  0.01e-3;        % modulation period in s
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
a       = 4*a0;      % Scattering length
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


[J,U,alpha,beta,gamma,epsilon,delta,chi] = OL2TB_parameters(nAtoms , s , lLat , omgPerp, omgLong , mCs , a )


%pi= 4.0*atan(1.0);
%L = 65;            % number of lattice sites
%N = 10e4;          % Number of particles, for normalization of the wavefunction

%Hubbard model parameters
%J = -1.0;
%U =  0.5;

%Harmonic trap: local energy shifts
% deltaE_i = epsilon*x_i^2+ delta
%epsilon = 100.0;
%delta   = 0.0;

%Dissipation
%gamma = 16.0

%%

%%%******************************************* 
% ode45 solver settings
%******************************************* 
options  = odeset('RelTol',3e-7,'AbsTol',1e-6,'Stats','off','MaxStep',1e-3);


psi      = zeros(2,1);
N_time   = 16;
tau      = zeros(N_time+1,1);
%epsilons = transpose(linspace(40,100,N_time+1));

Trap_period = transpose(linspace(20e-3,120e-3,N_time+1));
K           = 3.26;  % Argument of the Bessel function 
gamma_damping =  0.00;


for r=2:N_time+1
    %epsilon=epsilons(r);
    
    omgLong = 2.0*pi/Trap_period(r);
    [J,U,alpha,beta,gamma,epsilon,delta,chi] = OL2TB_parameters(nAtoms , s , lLat , omgPerp, omgLong , mCs , a )
    
    J_eff = J*besselj(0,K);
    epsilon_ = mCs*(omgLong*lLat)^2/Er;
    
    psi(1) = 0.0;    % initial position
    psi(2) = 0.5;    % initial momentum 
    
    %tSteps = linspace(0,200,256);
    tSteps = linspace(0,3.5/epsilon_^0.6,256);

    %propagation of the Schrodinger equation. 
    [TT,PSI_out] = ode45(@(t,y) DGLFUNC(J_eff,epsilon_,gamma_damping,t,y,0),tSteps,psi,options);

    psi_t = PSI_out;                       % take the last time step
    figure(2)
    plot(TT,psi_t(:,2));%,TT,psi_t(:,2));
    hold on
    title('SC');
    xlabel('time (ms)');
    ylabel('quasimomentum)');
    drawnow;
    [c,i]=min(abs(psi_t(:,2)-1.0));
    tau(r)=TT(i);
end

figure(4)
plot(1e3*Trap_period(2:end),1e3*(2*pi*hbar/Er)*tau(2:end),'.-','LineWidth',1);
axis([20 120 20 80])
hold on
title('Semiclassical model');
xlabel('Trap period T (ms)');
ylabel('Decat time (ms)');
drawnow;

%% 
function dpsi = DGLFUNC(J,epsilon_,gamma_damping,t,psi,RealImag)
    dpsi    = zeros(2,1);
    dpsi(1) = 2*pi*J*sin(pi*psi(2));    
    dpsi(2) = -epsilon_*psi(1)- gamma_damping*dpsi(1);
end


