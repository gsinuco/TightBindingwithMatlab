clear
clc

global L N pi ;

pi = 4.0*atan(1.0);
L  = 128;            % number of lattice sites
N  = 10e4;          % Number of particles, for normalization of the wavefunction


x  = linspace(-64,64,L);
psi = zeros(L);
psi = (psi+1);

Translation=translation(); % Definition of the translation operator: T|i> = |i+1>

J    = 30.0; % Tunneling rate
H_NN = zeros(L,L);
H_NN = NN_tunneling(J); %Near-Neighbour tunneling matrix

epsilon = 0.0125e-3;
delta   = 0.0e-3;
H_U2    = zeros(L,L);
H_U2    = onsite_shift(epsilon,delta); % On-site energy shift due to the trapping potential
                                       % delta E_i = epsilon*i^2 + delta


[V,D] =eig(H_NN+H_U2);  % Energy spectrum in the trapped OL

%transform Translation: Translation operator in the basis of eigenstates
Translation = transpose(V)*Translation*V;

figure(1)
subplot(1,2,1)
plot(x,diag(D), '.-')
subplot(1,2,2)
plot(x,V(:,1),'.-',x,V(:,2),'.-',x,V(:,3),'.-')
%plot(x,V(:,17),'-')
V_new = zeros(L,L);


V_new(:,1) = V(:,1);
V_new(:,L) = V(:,L);
%if(mod(L,2)==1) 
if(abs(D(1,1)-D(2,2))>1e-6)
    i_start = 1;
else
%if(mod(L,2)==0) 
    i_start = 0;
end
U_2PW=fromStationary2PlaneWaves(i_start); % In the periodic case, the calculated eigenfunctions
                                          % are real, and correspond to
                                          % cos()and sin().
                                          % This function produces the transformition the
                                          % wavefunctions to plane waves.
                                          
%here we do the transformation to plane-waves by hand.
for i=i_start:L/2-1
    if(i_start==1)
        V_new(:,2*i)   = (V(:,2*i)+1i*V(:,2*i+1))/sqrt(2);
        V_new(:,2*i+1) = (V(:,2*i)-1i*V(:,2*i+1))/sqrt(2);
    end
    if(i_start==0)
        V_new(:,2*i+1)   = (V(:,2*i+1)+1i*V(:,2*i+2))/sqrt(2);
        V_new(:,2*i+2) = (V(:,2*i+1)-1i*V(:,2*i+2))/sqrt(2);
    end
end 

V_new2 = V*U_2PW;

figure(2)
subplot(1,2,1)
plot(x,imag(V_new(:,1)),'.-',x,imag(V_new(:,2)),'.-')%,x,imag(V_new(:,6)),'.-')
subplot(1,2,2)
plot(x,real(V_new(:,1)),'.-',x,real(V_new(:,2)),'.-')%,x,real(V_new(:,6)),'.-')


figure(3)
subplot(1,2,1)
plot(x,imag(V_new2(:,1)),'.-',x,imag(V_new2(:,2)),'.-')%,x,imag(V_new(:,6)),'.-')
subplot(1,2,2)
plot(x,real(V_new2(:,1)),'.-',x,real(V_new2(:,2)),'.-')%,x,real(V_new(:,6)),'.-')

L_2 = floor(L/2);
for i=1:L
    k(i) = angle(V_new2(L_2+1,i)/V_new2(L_2,i));
end

[out,idx] = sort(k);

E_band = diag(D);

figure(4)
plot(k(idx),E_band(idx),'.-')

figure(5)
plot(x,abs(V_new2(:,1)),x,abs(V_new2(:,L)),'.',x,abs(V(:,1)),x,abs(V(:,L)))

%%
U   =  1.0e-2; % Hubbard mean-field on-site interaction
H_U = zeros(L,L);
H_U = onsite(U,psi);

epsilon = 0.1;
delta   = 0.0e-3;
H_U2 = zeros(L,L);
H_U2 = onsite_shift(epsilon,delta); % ON-site energy shift due to the trap
                                    % delt E_i= epsilon i^2 + gamma
alpha =  1.0e-1;
beta  = -alpha;
gamma = 0.0;
H_NN2 = zeros(L,L);
H_NN2 = NN_tunneling_site(alpha,beta,gamma); %Site-dependent correction to the tunneling
                                             % DELTA_I = alpha*i^2 + beta*i
                                             % + gamma

chi  = 1.0;
H_NLNNT = zeros(L,L);
H_NLNNT = NL_tunneling(chi,psi); % density-dependent tunneling rate
                                 % chi .....

H = zeros(L,L);
H = H_NN + 1.0*H_U + 0.0*H_NN2 + + 1.0*H_U2 + 0.0*H_NLNNT ;

%initial state
[V,D] =eig(H);
%+H_U2);%+H_U);
psi = V(:,1);
H_U = zeros(L,L);
H_U = onsite(U,psi);
%[V,D] =eig(H_NN+H_U2+H_U);
%psi = V(:,1);
figure(2)
plot(x,abs(psi).^2,x,abs(V(:,2)).^2,x,abs(V(:,3)).^2)   

xlabel('Position'); 
ylabel('Density');
%%
% Since H_NLNNT depends on the wavefunction
% we have to evaluate it again.
chi  = 1.0;
H_NLNNT = zeros(L,L);
H_NLNNT = NL_tunneling(chi,psi);
%%%
%******************************************* 
% ode45 solver settings
%******************************************* 
options=odeset('RelTol',3e-7,'AbsTol',1e-6,'Stats','off','MaxStep',1e-3);

%******************************************* 
% Find the ground state
%******************************************* 
Steps  = 128;
tSteps = linspace(0.0,10.0,Steps);

psi_x_2 = abs(psi).^2;  % Density psi_x                        
delta_psi = 10;         % Quantifier of the change of the psi
counter = 0;            % iterations counter
dt_gs   = 1e-1;         % time-step to find ground state

while (delta_psi > 1e-6)
    counter = counter+1
    %Imaginary propagation of the Schrodinger equation. 
    [TT,PSI_out] = ode45(@(t,y) DGLFUNC(L,t,y,U,J,epsilon,delta,0),[0 dt_gs/2 dt_gs],psi,options);
    
    psi = PSI_out(3,:);                % take the last time step
    psi = psi/sqrt(sum(abs(psi).^2));  % normalize
    psi_x_3 = abs(psi).^2;             % Density psi_x                       
    if mod(counter,10)==0
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
        H =  H_NN + 1.0*H_U2 + H_U;
        dpsi = -H*psi;
    end
    if(RealImag == 1) % Real time evolution
        H = -H_NN + 1.0*H_U2 + H_U;
        dpsi = 1i*H*psi;
    end
end

