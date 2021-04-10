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

%%%******************************************* 
% ode45 solver settings
%******************************************* 
options  = odeset('RelTol',3e-7,'AbsTol',1e-6,'Stats','off','MaxStep',1e-3);

%%%******************************************* 
% Energy Spectrum OL plus harmonic trap
%******************************************* 



function   V = V_trap(x,omega_trap,D,L):
    
    omega_trap2 = omega_trap*omega_trap;
    alpha       = 2*omega_trap2* D*L/(D-L) ;
    beta        = omega_trap2*D*(2*L-D)/((D-L)**2);
    
    V = zeros_like(x);
    
    if (D<L):
        dV_D = 2*omega_trap2*D/2;
        V_D  = omega_trap2*D*D/2;   
        a    = dV_D/(2*(D-L));
        b    = -dV_D*L/(D-L);
        c    = V_D3 - a*D*D  - b*D; 
    
        to_the_left  = np.argwhere(x < -D);
        to_the_right = np.argwhere(x >  D);
        central_section     = np.argwhere(np.abs(x)<=D);
        V[to_the_left]      = a*x[to_the_left]**2-b*x[to_the_left]+c;
        V[to_the_right]     = a*x[to_the_right]**2+b*x[to_the_right]+c;
        V[central_section]  = omega_trap2*x[central_section]**2/2;
    else:
        V = omega_trap2*x;
    end
end
            
    
    
function [E,DIAGONAL] = Spectrum_trapped_OL(L,N_x,D,omega_trap,V_0,phase_x):
    %k     = (1/(L*N_x))*(np.linspace(0,L*N_x,L*N_x))
    
    dx    = 1.0/N_x;

    H     = np.zeros([N_x*L,N_x*L],dtype=np.complex);
    x     = np.linspace(-L/2, L/2, N_x*L);


    V_trap_UnitCell = V_trap(x,omega_trap,D,L/2);


    V = np.diag(2/dx**2+0.5*V_0 + 0.5*V_0*(np.cos(2.0*np.pi*x+phase_x))+V_trap_UnitCell)# + 0.5*omega_trap*(x**2));
    H = V + np.diag(-1/dx**2 * np.ones(x.shape[0]-1),k=1) + np.diag(-1/dx**2 * np.ones(x.shape[0]-1),k=-1);

    dim = H.shape[0];
    H[dim-1,0] = H[0,0] - 2/dx^2;
    H[0,dim-1] = H[0,0] - 2/dx^2;

    E = np.linalg.eigvalsh(H);
    DIAGONAL = np.diag(V)+2/dx^2;
    
end