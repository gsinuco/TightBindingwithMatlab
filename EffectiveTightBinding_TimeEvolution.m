function [time_Spam,q,psi_q_2_vect,Atoms_Edge] = EffectiveTightBinding_TimeEvolution(index_)
% figure(3)
% subplot(2,1,1)
% time_Spam = linspace(0,modPeriod*nSaved,nSaved);
% imagesc(1e3*time_Spam,q, psi_q_2_vect');
% 
% subplot(2,1,2)
% time_Spam = linspace(0,modPeriod*nSaved,nSaved);
% plot(1e3*time_Spam,Atoms_Edge);


%************************************************************************
% Simulate driven system with GPE 1D
%
%************************************************************************


%******************************************* 
% constants
%******************************************* 
c    = 299792458;       % Speed of light [m/s]
kB   = 1.3806503e-23;   % Boltzmann constant [J/K]
hbar = 1.054571817e-34;  % Plank constant [Js]
muB  = 9.27400949e-24;  % Bohr magneton [J/T]
%muB  = 9.27400949e-28; % Bohr magneton [J/G]
g    = 9.80553;         % Acceleration of gravity [m/s^2]
u    = 1.6605e-27;      % Atomic mass unit [kg]
a0   = 0.5291772083e-10;% Bohr's radius in [m]
mCs  = 133*u;           % Mass of 133Cs [kg]
lLat = 1064.49463e-9;   % Lattice wavelength [m]
kLat = 2*pi/lLat;       % Lattice wavevektor [1/m]
omgRec = hbar*kLat^2/2/mCs;% Recoil frequency [Hz]
G    = 2*kLat;             % Brilloin zone width [1/m]
Er   = hbar^2*kLat^2/2/mCs; % Recoil energy [J]
dLat = lLat/2;           % Lattice site distance [m]

gamma_ = 0.01; %Global damping

%******************************************* 
% Modulation and time steps
%******************************************* 

F0 = 1.0*mCs*g;

modPeriod =  1e-3; %2*hbar*kLat/F0;                         % modulation period in s
modOmega  =  2*pi/modPeriod;
nPeriods  = 50;                           % Number of periods for simulation
Steps = 32;                                % Saved steps per BO
dt = modPeriod;                            % Single time step dt
nSaved = 1 + nPeriods*Steps;               % number saved wave functions

%******************************************* 
% Physical system, general parameters 
%******************************************* 
N = 501;                          % Number of Wells, must be odd
s = 12;                           % Lattice Depth in Er
a = 1.0*4*a0;                     % Scattering length

N_time      = 10;
Trap_period = transpose(linspace(20e-3,120e-3,N_time));
periodLong  = Trap_period(index_);              % Trap period in seconds
omgLong     = 2*pi/periodLong;     % Longitudinal trap frequency
omgPerp     = 2*pi * 300;          % radial, perpendicular trap freq.

nAtoms = 1E4;                     % Number of atoms
expTime = 40e-3;                  % End expansion time

% length scales
sigLong = sqrt(hbar/mCs/omgLong); % Longitudinal oscillator length
sigPerp = sqrt(hbar/mCs/omgPerp); % Radial oscillator length

% assume that vertical beam is strong, large radial frequency
omgMean = (omgPerp*omgPerp*omgLong)^(1/3);    % Geometric averaged trap frequency
sigMean = sqrt(hbar/mCs/omgMean);           % Geometric averaged oscillator length
mu      = hbar*omgMean/2*((15*nAtoms*a/sigMean)^(2/5));      % Chemical Potential for TF-approx.
rTFLong = sqrt(2*mu/mCs/omgLong/omgLong);   % Longitudinal TF-radius
% XX check

omgLatt = 2*omgRec*sqrt(s);                       % Harmonic frequency in one well
sigLatt = sqrt(hbar/mCs/omgLatt);                 % Oscillator length in one well

MaxExpqr = hbar*kLat/mCs*expTime;                 % Expansion faktor


%******************************************* 
% Initialise vectors
%******************************************* 
dq = 2*pi/(dLat*(N-1));
x  = [-dLat*(N-1)/2:dLat:dLat*(N-1)/2];  % Real space vector
q = [-dq*(N-1)/2:dq:dq*(N-1)/2]/kLat;     % Quasi momentum vector
qexp = [-dq*(N-1):dq:dq*(N-1)]/kLat;      % Doubled Quasi momentum vector

% save parameters
PSI_save = zeros(nSaved,N);                          % Empty matrix for wave functions Psi(t,x)
PSI_save_q = zeros(nSaved,N);                        % Empty matrix for wave functions Psi(t,q)
PSI_save_Exp = zeros(nSaved,N*2-1);                  % Empty matrix for wave functions Psi(t,qexp)

PSI2_save = zeros(nSaved,N);                          % Empty matrix for wave functions Psi(t,x)^2
PSI2_save_q = zeros(nSaved,N);                        % Empty matrix for wave functions Psi(t,q)^2
PSI2_save_Exp = zeros(nSaved,N*2-1);                  % Empty matrix for wave functions Psi(t,qexp)^2

mean_q = zeros(1,nSaved);                            % Empty matrix for mean momentum of Psi(t,qexp)
mean2_q = zeros(1,nSaved);                           % Empty matrix for mean momentum of Psi(t,qexp)


%******************************************* 
% Potential energy
%******************************************* 
g3D  = 4*pi*hbar^2*a/mCs;                      % Nonlinear interaction term
Vext = (0.5*mCs*omgLong^2*x.^2)';                % Potential energie  


%******************************************* 
% Initialise wave function (Gauss or TF)
%******************************************* 
psi_x = zeros(1,N);                                 % Start Wave function
PSI_Exp_Start = zeros(1,N);                         % Scaled Start Wave function

if (a == 0)  % Gauss
	psi_x = exp(-(x.^2)/2/(sigLong)^2);      
    psi_x = psi_x/sqrt(sum(abs(psi_x).^2));         % Normalize
    
    %PSI_Exp_Start = exp(-(q.^2)/2/(sigLong/MaxExpqr)^2);
    %PSI_Exp_Start = abs(PSI_Exp_Start/sqrt(sum(abs(PSI_Exp_Start).^2))).^2;    
    %PSI_Exp_Start = abs(PSI_Exp_Start).^2 / sum(abs(PSI_Exp_Start).^2);
    % XXX
    
else % TF profile
    
    % init
    psi_x = x*0;
    ind = find( (x>-rTFLong) & (x<rTFLong) ); % inside
    psi_x(ind) = sqrt(rTFLong^2-(x(ind)).^2);

    %PSI_Exp_Start = 0*x;
    %ind = find( (q>-rTFLong/MaxExpqr) & (q<rTFLong/MaxExpqr) ); % inside
    %PSI_Exp_Start(ind) = sqrt((rTFLong/MaxExpqr)^2-q(ind).^2);

    % normalize
    psi_x = psi_x/sqrt(sum(abs(psi_x).^2));         % Normalize
    %PSI_Exp_Start = abs(PSI_Exp_Start/sqrt(sum(abs(PSI_Exp_Start).^2))).^2;
    % XXX
end


%******************************************* 
% Calculate bands
%******************************************* 
nCoeff   = 40;                      % Number of coeffs >0
mCoeff   = -nCoeff:1:nCoeff;        % List of coeffs
E0       = zeros(1,N);              % First bloch band energies
for index_q = 1:N
    H = make_HMatrix(q(index_q));
    [cn, En] = eig(H);
    E0(index_q) = En(1);
end
K = abs(E0(1)-E0((N-1)/2+1))/4*Er;
%%
%******************************************* 
% ode45 solver settings
%******************************************* 
%options=odeset('RelTol',3e-8,'AbsTol',1e-9,'Stats','off','MaxStep',1e-8);
options=odeset('RelTol',3e-7,'AbsTol',1e-6,'Stats','off','MaxStep',1e-6);


%******************************************* 
% Ground state, imagingary time propagation
%******************************************* 

psi_x_2 = abs(psi_x).^2;      % Density psi_x                        
delta = 10;
counter = 0;
dt_gs   = 1e-4;                  % dt for finding ground state

B = g3D/((2*pi)^(3/2))/sigPerp^2/sigLatt; 

while (delta > 1e-6)
    counter = counter+1;
    [TT,PSI_out] = ode45(@DGLFUNC2,[0 dt_gs/2 dt_gs],psi_x,options);
    
    psi_x = PSI_out(3,:);                    % take the last time step
    psi_x = psi_x/sqrt(sum(abs(psi_x).^2));  % normalize
    psi_x_3 = abs(psi_x).^2;                 % Density psi_x                       
    if mod(counter,10)==0
        plot(x/dLat,psi_x_3);
        title('Ground state');
        xlabel('Position');
        ylabel('psi(x)');
        drawnow;
    end
    delta = sum(abs(psi_x_3-psi_x_2));              % Has the Wave function changed?
    psi_x_2 = psi_x_3;
end
 

    function DPSI2 = DGLFUNC2(t,PSI)
        UDD = -4*pi*hbar^2*0.5*a0/mCs/((2*pi)^(3/2))/sigPerp^2/sigLatt;

        psi_j   = PSI;
        n       = nAtoms*abs(psi_j).^2;
        
        psi_jm1 = [0; PSI(1:N-1); ];        
        n_jm1   = nAtoms*abs(psi_jm1).^2;
        
        psi_jp1 = [PSI(2:N); 0;];
        n_jp1   = nAtoms*abs(psi_jp1).^2;
        
        %DNLS    = sqrt(n) * sqrt( mCs*omgPerp^2*g3D/(sqrt(2*pi)*pi*sigLatt) );
        DNLS   = B*n;%g3D * n/((2*pi)^(3/2))/sigPerp^2/sigLatt;
        %DPNSE  = min(DNL,DNLS);        
        DPSI2=-1/hbar*(Vext+DNLS).*psi_j+K/hbar*(psi_jm1+psi_jp1)+...
              -1/hbar*UDD*psi_j.*(n_jm1+n_jp1);
        %DPSI2  = -1/hbar*(Vext+DNL).*psi_j+K/hbar*(psi_jm1+psi_jp1);        
    end


%******************************************* 
% prepare loop
%******************************************* 

width_q = zeros(1,nSaved);   % Empty matrix for momentum width of Psi(t,qexp)
%width2_q = zeros(nSaved);   % Empty matrix for momentum width of Psi(t,qexp)
 
psi_x = sqrt(psi_x_3);    
psi_x_2 = psi_x_3;
   
psi_q = fftshift(fft(psi_x));
psi_q = psi_q/sqrt(sum(abs(psi_q).^2));
psi_q_2 = abs(psi_q).^2;
 

figure(2)
% density in space
subplot(1,2,1);
plot(x/dLat,psi_x_2, x/dLat,psi_x_2,'r');
axis([-50 50 0 inf]);

subplot(1,2,2);
plot(q,psi_q_2);
title('Momentum in q');
xlabel('Momentum q (kLat)');

drawnow;


%psi_exp_2 = conv(PSI_Exp_Start,psi_q_2);
%psi_exp_2 = psi_exp_2/sum(psi_exp_2);
 
 
%PSI_save(1,:) = psi_x_2;
%PSI_save_q(1,:) = psi_q_2;
%PSI_save_Exp(1,:) = psi_exp_2;
 
%Ewx = sum(psi_exp_2.*qexp);
%Erw2x = sum(psi_exp_2.*((qexp-Ewx).^2));
%mean_q(1) = Ewx;
%width_q(1) = 2*sqrt(Erw2x);


% count number of "occupied" lattice sites
test=sort(psi_x_2,'descend');
latticesites=1;
while (sum(test(1:latticesites))<0.999)
    latticesites = latticesites+1;
end
fprintf('Number of occupied lattice sites %i.\n',latticesites);
 

%******************************************* 
% % Plot simulation start
% %******************************************* 
% subplot(3,2,2);
% plot(x/dLat,psi_x_2, x/dLat,PSI_save(1,:),'r');
% title(['Prop. time = ' num2str(0,4) ' ms / ' ...
%  num2str(0,3) ' %']);
% xlabel('Lattice site n');
% 
% subplot(3,2,4);
% plot(q,psi_q_2);
% title('momentum in q');
% xlabel('q/kLat');
% subplot(3,2,6);
% %plot(qexp,psi_exp_2);
% subplot(3,2,[1 3 5]); plot(0,width_q(1),'o');
% xlabel('Propagation time [ms]'); ylabel('moment2_k');
% drawnow;

% g3d = 4*pi*hbar^2*a/m;                  % Nonlinear interaction term
% %Vexto = (x*Fext+0.5*m*omega_long^2*x.^2)';   % Potential due to Gravity and harmonic trap


%******************************************* 
% loop over time
%******************************************* 

%g3D = -g3D/100;
    
A = sqrt(mCs*omgPerp^2*g3D/(sqrt(2*pi)*pi*sigLatt));;
UDD = -4*pi*hbar^2*0.5*a0/mCs/((2*pi)^(3/2))/sigPerp^2/sigLatt;

    function DPSI= DGLFUNC(t, PSI, tau)

        %Fmod = F0 * cos( modOmega * (tau + t) );
        %Fmod = min(F0,F0/(30*modPeriod)*(tau+t)) * cos( modOmega * (tau + t) );
        Vext = (0.5*mCs*omgLong^2*(x).^2)';%(x*Fmod+ 0.5*mCs*omgLong^2*(x).^2)';
        %Vext = (x*F0)';

        psi_j = PSI;
        n     = nAtoms*abs(psi_j).^2;
        
        psi_jm1=[0; PSI(1:N-1); ];
        n_jm1  = nAtoms*abs(psi_jm1).^2;
        
        psi_jp1=[PSI(2:N); 0;];
        n_jp1 = nAtoms*abs(psi_jp1).^2;
        
        %DNLS=sqrt(n)*sqrt(mCs*omgPerp^2*g3D/(sqrt(2*pi)*pi*sigLatt));        
        DNLS   = B*n;%g3D * n/((2*pi)^(3/2))/sigPerp^2/sigLatt;
        %DPNSE=min(DNL,DNLS);
        %DPSI=-i/hbar*(Vext+DNL).*psi_j+i/hbar*K*(psi_jm1+psi_jp1);       
        %c = sqrt(1-sigLatt^2/sigPerp^2);
        
        %UDD = -4*pi*hbar^2*0.5*a0/mCs/((2*pi)^(3/2))/sigP erp^2/sigLatt;
        
        DPSI=(-1i-gamma_)/hbar*(Vext+DNLS).*psi_j+1.0*(1i-gamma_)/hbar*K*besselj(0,3.26)*(psi_jm1+psi_jp1)+ ... 
             (-1i-gamma_)/hbar*UDD*psi_j.*(n_jm1+n_jp1);
    end


tSteps = linspace(0,modPeriod,Steps);
% tVect  = (0:nSaved)*modPeriod * 1e3;
% tStepsFull = linspace(0,modPeriod,2*Steps);
% 
% % buffer
% psi_x_vect   = zeros( 2*length(tSteps)-1, length(x) );
% psi_x_2_vect = psi_x_vect;
% psi_q_vect   = psi_x_vect;
% psi_q_2_vect = psi_x_vect;
% psi_exp_2_vect = zeros( 2*length(tSteps)-1, length(psi_exp_2) );
 
nSaved = 240;
psi_x_2_vect = zeros(nSaved,N); 
psi_q_2_vect = zeros(nSaved,N); 
Atoms_Edge   = zeros(1,nSaved);
time_Spam = linspace(0,modPeriod*nSaved,nSaved);

for counter = 1:nSaved

    T = modPeriod*(counter-1);                       % start time of modulation period
    
    [TT,PSI_out1] = ode45(@(t,y) DGLFUNC(t,y,T),tSteps,psi_x,options);
    
    psi_x = PSI_out1(end,:);
    psi_x = psi_x/sqrt(sum(abs(psi_x).^2));
    psi_x_2 = abs(psi_x).^2;                          % Density psi_x                       
    psi_q = fftshift(fft(psi_x));
    psi_q = psi_q/sqrt(sum(abs(psi_q).^2));
    psi_q_2 = abs(psi_q).^2;
    
    psi_x_2_vect(counter,:) = psi_x_2; 
    psi_q_2_vect(counter,:) = psi_q_2;     
    
    Atoms_Edge(counter) = sum(psi_q_2(1:(N-1)/4))+sum(psi_q_2(3*(N-1)/4:end));
    
%     if(mod(counter,10)==0)
%         figure(2)
%      % density in space
%         subplot(1,2,1);
%         plot(x/dLat,psi_x_2, x/dLat,psi_x_2,'r');
%         axis([-50 50 0 inf]);
%         title(sprintf('Time = %i ms / %.2fperc',T*1e3, counter/nSaved*100 ) );
%         
%         subplot(1,2,2);
%         plot(q,psi_q_2);
%         title('Momentum in q');
%         xlabel('Momentum q (kLat)');
% 
%         drawnow;
% 
%     end

%     psi_exp_2 = conv(PSI_Exp_Start,psi_q_2);
%     psi_exp_2 = psi_exp_2/sum(psi_exp_2);

    
    % get all images
%     for j1=1:length(tSteps)
%         psi_x_vect(j1,:) =  PSI_out1(j1,:);
%         psi_x_vect(j1,:) =  psi_x_vect(j1,:)/sqrt(sum(abs(psi_x_vect(j1,:)).^2));
%         psi_x_2_vect(j1,:) = abs( psi_x_vect(j1,:) ).^2;  
%         
%         psi_q_vect(j1,:) = fftshift(fft(psi_x_vect(j1,:)));
%         psi_q_vect(j1,:) = psi_q_vect(j1,:)/sqrt(sum(abs(psi_q_vect(j1,:)).^2));
%         psi_q_2_vect(j1,:) = abs(psi_q_vect(j1,:)).^2;        
%     
%         psi_exp_2_vect(j1,:) = conv(PSI_Exp_Start,psi_q_2_vect(j1,:))';
%         psi_exp_2_vect(j1,:) = psi_exp_2_vect(j1,:)/sum(psi_exp_2_vect(j1,:));
%     end    
%     
%     %---------- save step ---------
%     PSI2_save(counter,:) = psi_x_2;
%     PSI2_save_q(counter,:) = psi_q_2;
%     PSI2_save_Exp(counter,:) = psi_exp_2;
%     Ewx = sum(psi_exp_2.*qexp);
%     Erw2x = sum(psi_exp_2.*((qexp-Ewx).^2));
%     mean2_q(counter) = Ewx;
%     width2_q(counter) = 2*sqrt(Erw2x);

%  
%     [TT,PSI_out2] = ode45(@(t,y) DGLFUNC(t,y,T+modPeriod/2),tSteps,psi_x,options);
%     
%     psi_x = PSI_out2(end,:);
%     psi_x = psi_x/sqrt(sum(abs(psi_x).^2));
%     psi_x_2 = abs(psi_x).^2;                          % Density psi_x    
%     
%     psi_q = fftshift(fft(psi_x));
%     psi_q = psi_q/sqrt(sum(abs(psi_q).^2));
%     psi_q_2 = abs(psi_q).^2;
%     
%     psi_exp_2 = conv(PSI_Exp_Start,psi_q_2);
%     psi_exp_2 = psi_exp_2/sum(psi_exp_2);
% 
%     % get all images
%     for j1=1:Steps-1
%         psi_x_vect(j1+Steps,:) =  PSI_out2(j1+1,:);
%         psi_x_vect(j1+Steps,:) =  psi_x_vect(j1+Steps,:)/sqrt(sum(abs(psi_x_vect(j1+Steps,:)).^2));
%         psi_x_2_vect(j1+Steps,:) = abs( psi_x_vect(j1+Steps,:) ).^2;  
%         
%         psi_q_vect(j1+Steps,:) = fftshift(fft(psi_x_vect(j1+Steps,:)));
%         psi_q_vect(j1+Steps,:) = psi_q_vect(j1+Steps,:)/sqrt(sum(abs(psi_q_vect(j1+Steps,:)).^2));
%         psi_q_2_vect(j1+Steps,:) = abs(psi_q_vect(j1+Steps,:)).^2;  
%         
%         psi_exp_2_vect(j1+Steps,:) = conv(PSI_Exp_Start,psi_q_2_vect(j1+Steps,:));
%         psi_exp_2_vect(j1+Steps,:) = psi_exp_2_vect(j1+Steps,:)/sum(psi_exp_2_vect(j1+Steps,:));
%     end
    
%     %***************************
%     % save steps
%     %***************************
% 
%     PSI_save(counter,:) = psi_x_2;
%     PSI_save_q(counter,:) = psi_q_2;
%     PSI_save_Exp(counter,:) = psi_exp_2;   
%     
%     Ewx = sum(psi_exp_2.*qexp);
%     Erw2x = sum(psi_exp_2.*((qexp-Ewx).^2));
%     
%     mean_q(counter) = Ewx;
%     width_q(counter) = 2*sqrt(Erw2x);
%     
%     %***************************
%     % plot steps
%     %***************************
%     figure(2)
%     subplot(3,2,1); 
%     plot(tVect(1:counter), width_q(1:counter),'o');
%     xlabel('Propagation time (ms)'); ylabel('Momentum width (kLat)');
% 
%     % density in space
%     subplot(3,2,2);
%     plot(x/dLat,psi_x_2, x/dLat,PSI_save(1,:),'r');
%     axis([-50 50 0 inf]);
%     title(sprintf('Time = %i ms / %.2fperc',T*1e3, counter/nSaved*100 ) );
%     xlabel('Lattice site #');
%         
%     subplot(3,2,3);
%     imagesc(tStepsFull,q, psi_q_2_vect');
%     caxis([0 0.0120]);
%     xlabel('Propagation time (ms)');
%     %plot(qexp,psi_exp_2);
%     title('Expansion density');
%     %xlabel('Momentum qexp (kLat)');
% 
%     subplot(3,2,4);
%     plot(q,psi_q_2);
%     title('Momentum in q');
%     xlabel('Momentum q (kLat)');
%     
%     subplot(3,2,[5 6]);
%     imagesc(tVect(1:counter), qexp, PSI_save_Exp(1:counter,:)');
%     caxis([0 0.0120]);
%     xlabel('Propagation time (ms)');
%     ylim([-1 1]);
% 
%     
% %     subplot(3,2,5);
% %     qdensity2_image = 8e3*PSI2_save_Exp(1:plotCount,:)';
% %     qdensity2_image(qdensity_image>64) = 64;
% %     image(tVect*1e3,qexp,qdensity2_image);
% %     xlabel('Propagation time [ms]');
% 
%     drawnow;

end
%%


dq = 2*pi/(dLat*(N-1));
x  = [-dLat*(N-1)/2:dLat:dLat*(N-1)/2];  % Real space vector
q = [-dq*(N-1)/2:dq:dq*(N-1)/2]/kLat;     % Quasi momentum vector
qexp = [-dq*(N-1):dq:dq*(N-1)]/kLat;      % Doubled Quasi momentum vector

% figure(3)
% subplot(2,1,1)

% imagesc(1e3*time_Spam,q, psi_q_2_vect');
% 
% subplot(2,1,2)
% time_Spam = linspace(0,modPeriod*nSaved,nSaved);
% plot(1e3*time_Spam,Atoms_Edge);

%%
function H = make_HMatrix(q)
    h1    = (2*mCoeff+q).^2+ s/2;
    hDiag = diag(h1,0);
    h2    = ones(2*nCoeff,1)*(+s)/4;
    hOffs = diag(h2,1)+diag(h2,-1);
    H = hDiag + hOffs;
end 

end
