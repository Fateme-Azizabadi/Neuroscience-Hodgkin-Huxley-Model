
%%
clc;
clear;
%%
%HodgkinHuxley
Vrest = -60; % mV
dt = 0.01; % ms
Duration=20000;
T=ceil(Duration/dt);
totalTime = 150; % ms
C = 1; % uF/cm^2
% constants; values based on Table 1
E_Na = 55; % mV
E_K = -72; %mV
E_Leak = -49.4; % mV
g_Na = 120; % mS/cm^2
g_K = 36; % mS/cm^2
g_Leak = 0.3; % mS/cm^2
% Vector of timesteps
t = [0:dt:totalTime];

% Current input __ change this to see how different inputs affect the neuron
I_current = ones(1,length(t))*0.0;
I_current(50/dt:end) =7; % Input of 3 microA/cm2 beginning at 50 ms and steady until end of time period.
%I_current(50/dt:51.1/dt) = 6.3;%Part 4
%I_current(50/dt:end) = linspace(2,200,length(I_current(50/dt:end)));%Part7
%Part 9_Current
% Triangular Current
%I_=7;
%x =I_ .* sawtooth(2*pi*(1/50)*t,1/2)+ I_;
%I_current(50/dt:end) =x(50/dt:end);
% Pulse 
%I_current(50/dt:100/dt) =7;
%I_current(140/dt:end) =7;
%Chrip
%A=10;
%omega=0.02*pi;
%s_t=t.^2/4;
%y_t=A*cos(omega*t+s_t);
%I_current(50/dt:end) =y_t(1:length(I_current(50/dt:end)));
% Sin()
%A=10;
%omega=2*pi;
%y_t=A*sin(omega*t);
%I_current(50/dt:end) =y_t(1:length(I_current(50/dt:end)));

% initializing values
V(1) = Vrest; % membrane potential is starting at its resting state
% separate functions to get the alpha and beta values
[alphaM, betaM] = m_equations(V(1), Vrest);
[alphaN, betaN] = n_equations(V(1), Vrest);
[alphaH, betaH] = h_equations(V(1), Vrest);
% initializing gating variables to the asymptotic values when membrane potential
% is set to the membrane resting value based on equation 13
m(1) = (alphaM / (alphaM + betaM));
n(1) = (alphaN / (alphaN + betaN));
h(1) = (alphaH / (alphaH + betaH));
% repeat for time determined in totalTime , by each dt
for i = 1:length(t)
    % calculate new alpha and beta based on last known membrane potenatial
    [alphaN, betaN] = n_equations(V(i), Vrest);
    [alphaM, betaM] = m_equations(V(i), Vrest);
    [alphaH, betaH] = h_equations(V(i), Vrest);
    % conductance variables-computed separately to show how this
    % changes with membrane potential in one of the graphs
    conductance_K(i) = g_K*(n(i)^4);
    conductance_Na(i)=g_Na*(m(i)^3)*h(i);
    % retrieving ionic currents
    I_Na(i) = conductance_Na(i)*(V(i)-E_Na);
    I_K(i) = conductance_K(i)*(V(i)-E_K);
    I_Leak(i) = g_Leak*(V(i)-E_Leak);
    % Calculating the input
    Input = I_current(i) - (I_Na(i) + I_K(i) + I_Leak(i));
    % Calculating the new membrane potential
    V(i+1) = V(i) + Input* dt*(1/C);
    % getting new values for the gating variables
    m(i+1) = m(i) + (alphaM *(1-m(i)) - betaM * m(i))*dt;
    n(i+1) = n(i) + (alphaN *(1-n(i)) - betaN * n(i))*dt;
    h(i+1) = h(i) + (alphaH *(1-h(i)) - betaH * h(i))*dt;
end 
%%
%Part 1
% Special graph to show ionic current movement
Vrest = -60;
voltage = [-80:0.01:10];
for i = 1:length(voltage)
    [alphaN, betaN] = n_equations(voltage(i), Vrest);
    [alphaM, betaM] = m_equations(voltage(i), Vrest);
    [alphaH, betaH] = h_equations(voltage(i), Vrest);
    taum(i) = 1/(alphaM+betaM);
    taun(i) = 1/(alphaN+betaN);
    tauh(i) = 1/(alphaH+betaH);
    xm(i) = alphaM/(alphaM+betaM);
    xn(i) = alphaN/(alphaN+betaN);
    xh(i) = alphaH/(alphaH+betaH);
    aN(i) = alphaN;
    bN(i) = betaN;
    aM(i) = alphaM;
    bM(i) = betaM;
    aH(i) = alphaH;
    bH(i) = betaH;
end
% Steady State Value of m,n,h
figure('Name', ' Steady State Value');
plot(voltage, xh, voltage, xm, voltage, xn, 'LineWidth', 2);
legend('h_\infty','m_\infty','n_\infty');
title('Steady State Value');
xlabel('mV');
ylabel('x_\infty(V)');
grid on
% Time Constant
figure('Name', ' Time Constant');
plot(voltage,tauh,voltage,taum,voltage,taun,'LineWidth', 2);
legend('\tau_h','\tau_m','\tau_n');
title('Time Constant');
xlabel('mV');
ylabel('\tau_i(V)');
grid on
%%
% Part 2
%Action Potential
figure('Name', 'Membrane Potential')
plot(t(45/dt:end),V(45/dt:end-1), 'LineWidth', 2)
xlabel('Time (ms)')
ylabel('Voltage (mV)')
title('Action Potential')
grid on

%Action Potential
figure('Name', 'Membrane Potential')
plot(t(45/dt:80/dt),V(45/dt:80/dt), 'LineWidth', 2)
xlabel('Time (ms)')
ylabel('Voltage (mV)')
title('Action Potential')
grid on

%%
% Part 3
figure('Name', 'Membrane Potential vs input')
subplot(2,1,1)
plot(t(45/dt:end),V(45/dt:end-1), 'LineWidth', 2)
ylim([-70 50])
xlabel('Time (ms)')
ylabel('Voltage (mV)')
title('Action Potential')
grid on
subplot(2,1,2)
plot(t(45/dt:end),I_current(45/dt:end), 'r', 'LineWidth', 2)
xlabel('Time (ms)')
ylabel('Current (\muA)')
title('Input')
grid on 

%%
%Part 4
figure('Name', 'Membrane Potential vs input')
subplot(2,1,1)
plot(t(45/dt:80/dt),V(45/dt:80/dt), 'LineWidth', 2)
ylim([-70 50])
xlabel('Time (ms)')
ylabel('Voltage (mV)')
title('Action Potential')
grid on
subplot(2,1,2)
plot(t(45/dt:80/dt),I_current(45/dt:80/dt), 'r', 'LineWidth', 2)
xlabel('Time (ms)')
ylabel('Current (\muA)')
title('Input')
grid on

%%
% Part 5
% Calculation of Amplitude and frequency
[p,l]=findpeaks(V(45/dt:end-1),t(45/dt:end));
Amp=mean(p);% Calculation of Amplitude
diff_T=diff(l);
T_Voltage=mean(diff_T);
Freq=(1/T_Voltage)*1000;% Frequenc(Hz)

%Plotting
figure('Name', 'Membrane Potential vs input')
subplot(3,1,1)
plot(t(45/dt:end),V(45/dt:end-1),l,p,'o', 'LineWidth', 2)
ylim([-70 50])
xlabel('Time (ms)')
ylabel('Voltage (mV)')
title('Action Potential')
grid on
subplot(3,1,2)
plot(t(45/dt:end),I_current(45/dt:end), 'r', 'LineWidth', 2)
xlabel('Time (ms)')
ylabel('Current (\muA)')
title('Input')
grid on
ax=subplot(3,1,3);
srt={'current=7(\muA)',"Amplitude="+Amp+" (mV)","Frequency="+Freq+" (Hz)"};
text(0, 0.5,srt)
set (ax,'visible', 'off')
%text(60/dt, 4, )

%%
% Part 6

figure('Name', 'Membrane Potential vs input vs v-n')
subplot(3,1,1)
plot(t(45/dt:end),V(45/dt:end-1) ,'LineWidth', 2)
hold on
vt=[0:dt:totalTime];
vt(45/dt:end)=-55;
plot(t(45/dt:end),vt(45/dt:end),'--',color='m')
ylim([-70 50])
xlabel('Time (ms)')
ylabel('Voltage (mV)')
title('Action Potential')
grid on
subplot(3,1,2)
plot(t(45/dt:end),I_current(45/dt:end), 'r', 'LineWidth', 2)
xlabel('Time (ms)')
ylabel('Current (\muA)')
title('Input')
grid on

subplot(3,1,3)
plot(V(45/dt:end),n(45/dt:end), 'b', 'LineWidth', 2)
ylim([0 1])
xlabel('Membrance Potential (mV)')
ylabel('K^+ Activation Gate,n')
grid on

%%
% Part 8
figure('Name', 'Membrane Potential vs input vs v-n vs v-m vs v-h')

subplot(3,1,1)
plot(V(45/dt:end),n(45/dt:end), 'b', 'LineWidth', 2)
ylim([0 1])
xlabel('Membrance Potential (mV)')
ylabel('K^+ Activation Gate,n')
grid on
subplot(3,1,2)
plot(V(45/dt:end),m(45/dt:end), 'LineWidth', 2)
ylim([0 1])
xlabel('Membrance Potential (mV)')
ylabel('Na^+ Activation Gate,m')
grid on
subplot(3,1,3)
plot(V(45/dt:end),h(45/dt:end), 'LineWidth', 2)
ylim([0 1])
xlabel('Membrance Potential (mV)')
ylabel('Na^+ Inactivation Gate,h')
grid on

%%
% Part 9
for j=1:10:300
    I_current(50/dt:end) = j;%Part 4
    % initializing values
    V(1) = Vrest; % membrane potential is starting at its resting state
    % separate functions to get the alpha and beta values
    [alphaM, betaM] = m_equations(V(1), Vrest);
    [alphaN, betaN] = n_equations(V(1), Vrest);
    [alphaH, betaH] = h_equations(V(1), Vrest);
    % initializing gating variables to the asymptotic values when membrane potential
    % is set to the membrane resting value based on equation 13
    m(1) = (alphaM / (alphaM + betaM));
    n(1) = (alphaN / (alphaN + betaN));
    h(1) = (alphaH / (alphaH + betaH));
    % repeat for time determined in totalTime , by each dt
    for i = 1:length(t)
        % calculate new alpha and beta based on last known membrane potenatial
        [alphaN, betaN] = n_equations(V(i), Vrest);
        [alphaM, betaM] = m_equations(V(i), Vrest);
        [alphaH, betaH] = h_equations(V(i), Vrest);
        % conductance variables-computed separately to show how this
        % changes with membrane potential in one of the graphs
        conductance_K(i) = g_K*(n(i)^4);
        conductance_Na(i)=g_Na*(m(i)^3)*h(i);
        % retrieving ionic currents
        I_Na(i) = conductance_Na(i)*(V(i)-E_Na);
        I_K(i) = conductance_K(i)*(V(i)-E_K);
        I_Leak(i) = g_Leak*(V(i)-E_Leak);
        % Calculating the input
        Input = I_current(i) - (I_Na(i) + I_K(i) + I_Leak(i));
        % Calculating the new membrane potential
        V(i+1) = V(i) + Input* dt*(1/C);
        % getting new values for the gating variables
        m(i+1) = m(i) + (alphaM *(1-m(i)) - betaM * m(i))*dt;
        n(i+1) = n(i) + (alphaN *(1-n(i)) - betaN * n(i))*dt;
        h(i+1) = h(i) + (alphaH *(1-h(i)) - betaH * h(i))*dt;
    end 
    [p,l]=findpeaks(V(45/dt:end-1),t(45/dt:end));
    Amp=mean(p);% Calculation of Amplitude
    diff_T=diff(l);
    T_Voltage=mean(diff_T);
    Freq=(1/T_Voltage)*1000;% Frequenc(Hz)
    F=ceil(Freq);% T=1 s
    if j>160
        F_Rate(floor(j/10)+1)=0;
    else
        F_Rate(floor(j/10)+1)=F;
    end
    
end
I=[1:10:300];
%Plotting
figure('Name', 'F-I Curve')
stem(I,F_Rate,'LineWidth', 2)
xlim([0 200])
xlabel('Current (\muA) ')
ylabel('Firing Rate (number of spikes in one second)')
title('F-I Curve')
grid on







