Atoms_Edge_TimeEvol = zeros(10,240);
for period_index = 1:10
    [time_Spam,q,psi_q_2_vect,Atoms_Edge] = EffectiveTightBinding_TimeEvolution(period_index);
    %[time_Spam,q,psi_q_2_vect,Atoms_Edge] = EffectiveTightBinding_TimeEvolution_V2(period_index);
    %psi_q_2_TimeEvol(period_index,:)    = psi_q_2_vect;
    Atoms_Edge_TimeEvol(period_index,:) = Atoms_Edge;
end
%%

tau     = zeros(1,10);
periods = zeros(1,10);
Trap_period = transpose(linspace(20e-3,120e-3,10));
for period_index = 1:10
    periods(period_index) = Trap_period(period_index);              % Trap period in seconds
    counter = 1;
    found = 0;
    while found==0
        
        if Atoms_Edge_TimeEvol(period_index,counter) > 0.5
        %if 1-Atoms_Edge_TimeEvol(period_index,counter) > 1-Atoms_Edge_TimeEvol(period_index,counter-1)
            found = 1;
            tau(period_index) = counter
        else
            found = 0;
        end
        counter = counter + 1;
    end
end


%%
Atoms_Edge_TimeEvol_FT = zeros(10,160);
for period_index = 1:10
    [time_Spam_FT,q,psi_q_2_vect,Atoms_Edge_FT] = EffectiveTightBinding_FiniteTemp_TimeEvolution_V2(period_index);
    %psi_q_2_TimeEvol(period_index,:)    = psi_q_2_vect;
    Atoms_Edge_TimeEvol_FT(period_index,:) = Atoms_Edge_FT;
end

%%
tau_FT     = zeros(1,10);
periods_FT = zeros(1,10);

for period_index = 1:10
    periods_FT(period_index) = 120e-3/period_index;              % Trap period in seconds
    counter = 1;
    found = 0;
    while found==0
        
        if Atoms_Edge_TimeEvol_FT(period_index,counter) > 0.5
        %if 1-Atoms_Edge_TimeEvol(period_index,counter) > 1-Atoms_Edge_TimeEvol(period_index,counter-1)
            found = 1;
            tau_FT(period_index) = counter
        else
            found = 0;
        end
        counter = counter + 1;
    end
end

%%
timeEvol_SC = zeros(10,4096);
[periods_SC,timeEvol_SC,tau_ms] = SC_V3(1.0,0.0);

timeEvol_SC = zeros(10,4096);
[periods_SC,timeEvol_SC_gamma,tau_ms_gamma] = SC_V3(1.0,0.55);

%%
figure(3)
subplot(2,1,1)
imagesc(1e3*time_Spam_FT,q, psi_q_2_vect');

figure(4)
subplot(10,1,1)
plot(1e3*time_Spam_FT,Atoms_Edge_TimeEvol_FT(1,:),'LineWidth',2);
legend('T=120 ms')
ylabel('Fraction');
subplot(10,1,2)
plot(1e3*time_Spam_FT,Atoms_Edge_TimeEvol_FT(2,:),'LineWidth',2);
legend('T=108 ms')
ylabel('Fraction');
subplot(10,1,3)
plot(1e3*time_Spam_FT,Atoms_Edge_TimeEvol_FT(3,:),'LineWidth',2);
legend('T=96 ms')
ylabel('Fraction');
subplot(10,1,4)
plot(1e3*time_Spam_FT,Atoms_Edge_TimeEvol_FT(4,:),'LineWidth',2);
legend('T=84 ms')
ylabel('Fraction');
subplot(10,1,5)
plot(1e3*time_Spam_FT,Atoms_Edge_TimeEvol_FT(5,:),'LineWidth',2);
legend('T=72 ms')
ylabel('Fraction');
subplot(10,1,6)
plot(1e3*time_Spam_FT,Atoms_Edge_TimeEvol_FT(6,:),'LineWidth',2);
legend('T=60 ms')
ylabel('Fraction');
subplot(10,1,7)
plot(1e3*time_Spam_FT,Atoms_Edge_TimeEvol_FT(7,:),'LineWidth',2);
legend('T=48 ms')
ylabel('Fraction');
subplot(10,1,8)
plot(1e3*time_Spam_FT,Atoms_Edge_TimeEvol_FT(8,:),'LineWidth',2);
legend('T=36 ms')
ylabel('Fraction');
subplot(10,1,9)
plot(1e3*time_Spam_FT,Atoms_Edge_TimeEvol_FT(9,:),'LineWidth',2);
legend('T=24 ms')
ylabel('Fraction');
subplot(10,1,10)
plot(1e3*time_Spam_FT,Atoms_Edge_TimeEvol_FT(10,:),'LineWidth',2);
legend('T=12 ms')
xlabel('Propagation time [ms]'); 
ylabel('Fraction');
%%
figure(3)
subplot(1,1,1)
imagesc(1e3*time_Spam,q, psi_q_2_vect');
xlabel('Propagation time [ms]'); 
ylabel('quasimomentum');
%%
figure(4)
subplot(10,1,1)
hold on
plot(1e3*time_Spam,Atoms_Edge_TimeEvol(1,:),'LineWidth',2);
r=1;
tSteps = linspace(0,240,4096);
plot(tSteps,-timeEvol_SC(r,:),'LineWidth',2);
legend('T=120 ms')
ylabel('Fraction');
subplot(10,1,2)
hold on
plot(1e3*time_Spam,Atoms_Edge_TimeEvol(2,:),'LineWidth',2);
plot(tSteps,-timeEvol_SC(2,:),'LineWidth',2);
legend('T=108 ms')
ylabel('Fraction');
subplot(10,1,3)
hold on
r =3;
plot(1e3*time_Spam,Atoms_Edge_TimeEvol(3,:),'LineWidth',2);
plot(tSteps,-timeEvol_SC(r,:),'LineWidth',2);
legend('T=96 ms')
ylabel('Fraction');
subplot(10,1,4)
hold on
r = 4;
plot(1e3*time_Spam,Atoms_Edge_TimeEvol(4,:),'LineWidth',2);
plot(tSteps,-timeEvol_SC(r,:),'LineWidth',2);
legend('T=84 ms')
ylabel('Fraction');
subplot(10,1,5)
hold on
r = 5;
plot(1e3*time_Spam,Atoms_Edge_TimeEvol(5,:),'LineWidth',2);
plot(tSteps,-timeEvol_SC(r,:),'LineWidth',2);
legend('T=72 ms')
ylabel('Fraction');
subplot(10,1,6)
hold on
r =6;
plot(1e3*time_Spam,Atoms_Edge_TimeEvol(6,:),'LineWidth',2);
plot(tSteps,-timeEvol_SC(r,:),'LineWidth',2);
legend('T=60 ms')
ylabel('Fraction');
subplot(10,1,7)
hold on 
r =7
plot(1e3*time_Spam,Atoms_Edge_TimeEvol(7,:),'LineWidth',2);
plot(tSteps,-timeEvol_SC(r,:),'LineWidth',2);
legend('T=48 ms')
ylabel('Fraction');
subplot(10,1,8)
hold on
r = 8;
plot(1e3*time_Spam,Atoms_Edge_TimeEvol(8,:),'LineWidth',2);
plot(tSteps,-timeEvol_SC(r,:),'LineWidth',2);
legend('T=36 ms')
ylabel('Fraction');
subplot(10,1,9)
hold on 
r = 9;
plot(1e3*time_Spam,Atoms_Edge_TimeEvol(9,:),'LineWidth',2);
plot(tSteps,-timeEvol_SC(r,:),'LineWidth',2);
legend('T=24 ms')
ylabel('Fraction');
subplot(10,1,10)
hold on 
r = 10;
plot(1e3*time_Spam,Atoms_Edge_TimeEvol(10,:),'LineWidth',2);
plot(tSteps,-timeEvol_SC(r,:),'LineWidth',2);
legend('T=12 ms')
xlabel('Propagation time [ms]'); 
ylabel('Fraction');

%%

data_Frequency = [8.2 11.9 16.9 20.6 23.9 26.9 29.2 31.6 33.8];
data_Period    = 1./data_Frequency;
data_Tau       = [1.473832e-02 1.656333e-02 2.059415e-02 2.335221e-02 2.715311e-02  ...
                  3.396861e-02 3.269406e-02 3.501393e-02 3.866137e-02];
%data_Period = [35 42 50 60 80 120];
%data_Tau    = [29 37 43 48 60 67];

%%

figure(5)
hold on
plot(1e3*data_Period,1./data_Tau,'.-','LineWidth',2,'MarkerSize',24);
plot(1e3*periods,tau, 'o-','LineWidth',2,'MarkerSize',8,'Color','g');
%plot(1e3*periods_FT,tau_FT-40, 'x-','LineWidth',2,'MarkerSize',8);
plot(periods_SC,tau_ms,'x-','LineWidth',2,'MarkerSize',8,'Color','b');
plot(periods_SC,tau_ms_gamma,'x-','LineWidth',2,'MarkerSize',8,'Color','r');
axis([20 120 20 80])
%title('Semiclassical Model');
xlabel('Trap period T (ms)');
ylabel('Decay time (ms)');
%legend('Experiment','Semiclassical','Semiclassical (damped)' )
legend('Experiment','1D-GPE','Semiclassical','Semiclassical (damped)')
drawnow;


%%

            
