%%
clear
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       PHASED ARRAY ANTENNA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('./Functions/');

% DATA
%data=;
%T_window=;

%MODULATION AND ENCODING
alpha_RS=255/223;
alpha_conv=6;
alpha_enc= alpha_conv*alpha_RS;
alpha_mod=1; 


TM_rate=6e3*alpha_mod/alpha_enc; %planning operations SARA
SD_rate=29000*alpha_mod/alpha_enc;   %104e3; %messenger site
R=TM_rate+SD_rate;

R_real=R*alpha_enc/alpha_mod;
%R=data/T_window;
f=8.4e9; %Hz
P_in=42; %tabelle sspa

mu_amp=0.323; %from graph
%BER=1e-5; %literature

P_tx=mu_amp*P_in;
P_tx=10*log10(P_tx);

% % phased array antenna
% wg_width=0.17/14.5 + 9*0.17/29; 
% wg_length=0.777;
% A_eff = wg_length*wg_width*8;

%GAIN
c=299792458;
lambda=c/f;
G_ph=10*log10(8) + 17; %https://www.antenna-theory.com/antennas/aperture/slottedwaveguide2.php

L_cable= -1; %dB

%GROUND
Drx=34;
mu_rx=0.55;
Grx=10*log10(pi^2*Drx^2*mu_rx/lambda^2);
theta_rx=65.3*lambda/Drx;

%Free space loss
r_e=astroConstants(2)*1e3;
r_m=58e9; %average value
r=1.553e11; %r_m+r_e;
L_space=20*log10(lambda/(4*pi*r));

%point losses
L_point=-0.1; %dB DSN datasheet
etarx=sqrt(L_point/(-12))*theta_rx;

%Atmospheric loss 
L_atm=-0.85; %dB from graph 

% Effective Isotropic Radiated Power
EIRP=P_tx+G_ph+L_cable

%Receiver power
P_rx = EIRP + Grx + L_space + L_point + L_atm

Kb=1.38e-23; %boltzman
Ts = 21; %Kelvin
N0 = 10*log10(Kb*Ts);

%Energy per bit to noise density
EbN0=P_rx-N0-10*log10(R_real)

% Modulation index
beta=deg2rad(68); %Mission oerview santo

%Carrier modulation index reduction
P_mod_loss=20*log10(cos(beta));

%carrier power
P_carrier= P_rx + P_mod_loss

SNR_carrier = 13; %dB 10 di minima +3

B = 10^(0.1*(- SNR_carrier + P_carrier - N0) )


%% 
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       MEDIUM GAIN ANTENNA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%MODULATION AND ENCODING
alpha_RS=255/223;
alpha_conv=6;
alpha_enc= alpha_conv*alpha_RS;
alpha_mod=1; 


% DATA
%data=;
%T_window=;
%TM_rate=; %planning operations SARA
%SD_rate=;   %104e3; %messenger site
R=3500*alpha_mod/alpha_enc;   %TM_rate+SD_rate;


R_real=R*alpha_enc/alpha_mod;
%R=data/T_window;
f=8.4e9; %Hz
P_in=42; %76.6; 

mu_amp=0.346; %from graph
BER=1e-5; %literature

P_tx=mu_amp*P_in;
P_tx=10*log10(P_tx);

%GAIN
c=299792458;
lambda=c/f;
G_ph=10*log10(0.8^2/2)+20.2; 

L_cable= -1; %dB

%GROUND
Drx=34;
mu_rx=0.55;
Grx=10*log10(pi^2*Drx^2*mu_rx/lambda^2);
theta_rx=65.3*lambda/Drx;

%Free space loss
r_e=astroConstants(2)*1e3;
r_m=58e9; %average value
r=1.553e11; %r_m+r_e;
L_space=20*log10(lambda/(4*pi*r));

%point losses
L_point=-0.1; %dB DSN datasheet
etarx=sqrt(L_point/(-12))*theta_rx;

%Atmospheric loss 
L_atm=-0.80; %dB from graph %internet sulle dsn della nasa BOB

% Effective Isotropic Radiated Power
EIRP=P_tx+G_ph+L_cable

%Receiver power
P_rx = EIRP + Grx + L_space + L_point + L_atm

Kb=1.38e-23; %boltzman
Ts = 21; %Kelvin
N0 = 10*log10(Kb*Ts);

%Energy per bit to noise density
EbN0=P_rx-N0-10*log10(R_real)

% Modulation index
beta=deg2rad(68); %Mission overview santo

%Carrier modulation index reduction
P_mod_loss=20*log10(cos(beta));

%carrier power
P_carrier= P_rx + P_mod_loss

SNR_carrier = 13; %dB 10 di minima +3

B = 10^(0.1*(- SNR_carrier + P_carrier - N0) )

%% 
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       LOW GAIN ANTENNA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MODULATION AND ENCODING
alpha_RS=255/223;
alpha_conv=6;
alpha_enc= alpha_conv*alpha_RS;
alpha_mod=1; 

% DATA
%data=;
%T_window=;
%TM_rate=; %planning operations SARA
%SD_rate=;   %104e3; %messenger site
R=350*alpha_mod/alpha_enc;   %TM_rate+SD_rate;

R_real=R*alpha_enc/alpha_mod;
f=8.4e9; %Hz
P_in=42; %76.6; 

mu_amp=0.346; %from graph
BER=1e-5; %literature

P_tx=mu_amp*P_in;
P_tx=10*log10(P_tx);

%GAIN
c=299792458;
lambda=c/f;
G_ph=6;
D=lambda*10^((G_ph-7)/20) ;

L_cable= -1; %dB

%GROUND
Drx=34;
mu_rx=0.55;
Grx=10*log10(pi^2*Drx^2*mu_rx/lambda^2);
theta_rx=65.3*lambda/Drx;

%Free space loss
r=1.553e11; %r_m+r_e;
L_space=20*log10(lambda/(4*pi*r));

%point losses
L_point=-0.1; %dB DSN datasheet
etarx=sqrt(L_point/(-12))*theta_rx;

%Atmospheric loss 
L_atm=-0.80; %dB from graph %internet sulle dsn della nasa BOB

% Effective Isotropic Radiated Power
EIRP=P_tx+G_ph+L_cable;

%Receiver power
P_rx = EIRP + Grx + L_space + L_point + L_atm;

Kb=1.38e-23; %boltzman
Ts = 221; %Kelvin
N0 = 10*log10(Kb*Ts);

%Energy per bit to noise density
EbN0=P_rx-N0-10*log10(R_real);

% Modulation index
beta=deg2rad(68); %Mission oerview santo

%Carrier modulation index reduction
P_mod_loss=20*log10(cos(beta));

%carrier power
P_carrier= P_rx + P_mod_loss;

SNR_carrier = 13; %dB 10 di minima +3

B = 10^(0.1*(- SNR_carrier + P_carrier - N0) );



