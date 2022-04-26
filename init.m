clear;              % Removes variables from the workspace
clc;                % Clear the screen
close all;          % Close the figures

CPU_clock = 200e6; % CPU clock frequency
Sys_clock_div = 2; % From hardware setting Auto set PLL
fsw = 10e3; % Switching frquency
 % Sampling period
PWM_timer_period = ...
    CPU_clock/Sys_clock_div/fsw/2; % Number 2 comes from Up-down carrier

%% Rate limiter values [rad/s]

max_rate = 100;
min_rate = -100;


%% Motor parameters (Teknic2310P)

%Rated values
U_N= 40;                    %V                      % Rated voltage
I_rated= 7.1;               %A                      % Rated current (phase-peak)
n_N  = 6000;                %rpm                    % Rated speed
p = 4;                      %                       % Number of pole pairs     
fN = 6000*p/60;             %Hz                     % Rated frequency

Ke= 4.64;                   %Bemf Const             % Vpk_LL/krpm
Kt= 0.0384;                 %Nm/A                   % Torque constant

psi_f = single((Ke)/(sqrt(3)*2*pi*1000*p/60));              % PM flux amplitude constant computed from Ke
psi_s = single(U_N/sqrt(3)*sqrt(2)/(2*pi*fN));              % Stator flux linkage
T_rated    = single((3/2)*p*psi_f*I_rated);                 % Get T_rated from I_rated

%Power supply Udc = 24V
U_m = 24/sqrt(3);           %V                      % Amplitude of phase voltage
n_n = U_m/psi_s*60/(2*pi)/p;%rpm                    % Maximum speed with constant magnetic flux linkage
w_Mn = U_m/psi_s/p;         %rads^-1                % Mechanical angular speed
w_mn = U_m/psi_s;
%Motor electrical and mechanical parameters
Rs = 0.36;                  %Ohm                    % Stator resistance  
Ld= 0.2e-3;                 %H                      % D-axis inductance value
Lq = 0.2e-3;                %H                      % Q-axis inductance value
Ls = 0.2e-3;                %H                      % Synchronous inductance of SPM(Ls=Ld=Lq)

J= 7.061551833333e-6;       %Kg-m2                  % Inertia in SI units
B= 2.636875217824e-6;       %Kg-m2/s                % Friction Co-efficient


% PositionOffset = 0.078;     %PU position            % Position Offset--> from old model
PositionOffset = 3.512;       % Measured position offset
QEPSlits       = 1000;      %                       % QEP Encoder Slits

%% Controller parameters
Ts = 1/(2*fsw);         % Sampling period, switching frequency fsw = 1/(2*Ts)
alphac = 2*pi*300;      % Current-controller bandwidth
alphas = 0.1*alphac;    % Speed-controller bandwidth
a = 2*pi*0.01;          % First order filter constant  