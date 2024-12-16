%% DATA RATE EBN0 LOW GAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   DATA RATE EBN0LOW ANTENNA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

addpath('./Functions/');

R_real=10:0.1:150;
EbN0=[];
for i=1:length(R_real)
    
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
    L_atm=-0.80; %dB from graph 
    
    % Effective Isotropic Radiated Power
    EIRP=P_tx+G_ph+L_cable;
    
    %Receiver power
    P_rx = EIRP + Grx + L_space + L_point + L_atm;
    
    Kb=1.38e-23; %boltzman
    Ts = 21; %Kelvin
    N0 = 10*log10(Kb*Ts);
    
    %Energy per bit to noise density
    EbN0_i=P_rx-N0-10*log10(R_real(i));

    EbN0=[EbN0; EbN0_i];
end

plot(R_real, EbN0, 'Color',[0.4660 0.6740 0.1880], 'LineWidth', 1)
xlabel('R_{real} [bps]')
ylabel('E_{b}/N_{0} [dB]')
grid on
hold on
yline(1.5, '--','Color',[0.4940 0.1840 0.5560], 'LineWidth', 1 )
hold on
yline(4.5, 'Color',[0.6350 0.0780 0.1840], 'LineWidth', 1)
hold on
xline(75.35, '--k', 'LineWidth', 1)
legend('E_{b}/N_{0}', 'E_{b}/N_{0}_{min}', 'E_{b}/N_{0}_{min} + margin', 'Maximum data rate allowed', 'Location', 'northeast' )
title('E_{b}/N_{0} in function of Data Rate for LGA')


%% DATA RATE EBN0 MEDIUM GAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   DATA RATE EBN0 MEDIUM ANTENNA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

R_real=[10:0.1:1000];
EbN0=[];
for i=1:length(R_real)
    
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
    G_ph=10*log10(0.8^2/2)+20.2; %https://www.antenna-theory.com/antennas/aperture/slottedwaveguide2.php
    
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
    L_atm=-0.80; %dB from graph 
    
    % Effective Isotropic Radiated Power
    EIRP=P_tx+G_ph+L_cable;
    
    %Receiver power
    P_rx = EIRP + Grx + L_space + L_point + L_atm;
    
    Kb=1.38e-23; %boltzman
    Ts = 21; %Kelvin
    N0 = 10*log10(Kb*Ts);
    
    %Energy per bit to noise density
    EbN0_i=P_rx-N0-10*log10(R_real(i));
    EbN0=[EbN0; EbN0_i];

end

plot(R_real, EbN0, 'Color',[0.4660 0.6740 0.1880], 'LineWidth', 1)
xlabel('R_{real} [bps]')
ylabel('E_{b}/N_{0} [dB]')
grid on
hold on
yline(1.5, '--','Color',[0.4940 0.1840 0.5560], 'LineWidth', 1 )
hold on
yline(4.5, 'Color',[0.6350 0.0780 0.1840], 'LineWidth', 1)
hold on
xline(634.212, '--k', 'LineWidth', 1)
legend('E_{b}/N_{0}', 'E_{b}/N_{0}_{min}', 'E_{b}/N_{0}_{min} + margin', ...
    'Maximum data rate allowed', 'Location', 'northeast' )
title('E_{b}/N_{0} in function of Data Rate for MGA')

%% DATA RATE EBN0 PHASED ARRAY ANTENNA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   DATA RATE EBN0 PHASED ARRAY ANTENNA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

R_real=[1000:1:10000];
EbN0=[];
for i=1:length(R_real)
      
    %R=data/T_window;
    f=8.4e9; %Hz
    P_in=42; %76.6; 

    mu_amp=0.323; %from graph
    BER=1e-5; %literature
    
    P_tx=mu_amp*P_in;
    P_tx=10*log10(P_tx);
    
    %GAIN
    c=299792458;
    lambda=c/f;
    G_ph=10*log10(8) + 17; %https://www.antenna-theory.com/antennas/aperture/slottedwaveguide2.php
    
    L_cable=-1; %dB
    
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
    L_atm=-0.80; %dB from graph 
    
    % Effective Isotropic Radiated Power
    EIRP=P_tx+G_ph+L_cable;
    
    %Receiver power
    P_rx = EIRP + Grx + L_space + L_point + L_atm;
    
    Kb=1.38e-23; %boltzman
    Ts = 21; %Kelvin
    N0 = 10*log10(Kb*Ts);
    
    %Energy per bit to noise density
    EbN0_i=P_rx-N0-10*log10(R_real(i));
        EbN0=[EbN0; EbN0_i];


end

plot(R_real, EbN0, 'Color',[0.4660 0.6740 0.1880], 'LineWidth', 1)
xlabel('R_{real} [bps]')
ylabel('E_{b}/N_{0} [dB]')
grid on
hold on
yline(1.5, '--','Color',[0.4940 0.1840 0.5560], 'LineWidth', 1 )
hold on
yline(4.5, 'Color',[0.6350 0.0780 0.1840], 'LineWidth', 1)
hold on
xline(7084.4, '--k', 'LineWidth', 1)
legend('E_{b}/N_{0}', 'E_{b}/N_{0}_{min}', 'E_{b}/N_{0}_{min} + margin', ...
    'Maximum data rate allowed', 'Location', 'northeast' )
title('E_{b}/N_{0} in function of Data Rate for PAA')

%% RAGGIO - LOW GAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   RAGGIO - LOW GAIN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%MODULATION AND ENCODING
alpha_RS=255/223;
alpha_conv=6;
alpha_enc= alpha_conv*alpha_RS;
alpha_mod=1; 


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

   %point losses
L_point=-0.1; %dB DSN datasheet
etarx=sqrt(L_point/(-12))*theta_rx;
    
    %Atmospheric loss 
    L_atm=-0.80; %dB from graph 
    
    % Effective Isotropic Radiated Power
    EIRP=P_tx+G_ph+L_cable;

    Kb=1.38e-23; %boltzman
    Ts = 21; %Kelvin
    N0 = 10*log10(Kb*Ts);

mu_S=astroConstants(4);

% Initial date
date_i=date2mjd2000([2011 03 18 0 0 0]);

% Final date
date_f=date2mjd2000([2015 04 30 0 0 0]);

date=date_i:1:date_f;

R_real=[];
r=[];

for i=1:length(date)
    

    % Mercury position
    [kep_m,ksun_m] = uplanet(date(i), 1);
    [r_m,v_m] = kep2car(kep_m(1), kep_m(2), kep_m(3), kep_m(4), kep_m(5), kep_m(6), mu_S);

    % Earth position
    [kep_e,ksun_e] = uplanet(date(i), 3);
    [r_e,v_e] = kep2car(kep_e(1), kep_e(2), kep_e(3), kep_e(4), kep_e(5), kep_e(6), mu_S);

    % Earth - Mercury distance
    r_i=r_m-r_e;
    r_i_norm=norm(r_i)*1e3;

    % Space loss
    L_space=20*log10(lambda/(4*pi*r_i_norm));

    % Received power
    P_rx_i = EIRP + Grx + L_space + L_point + L_atm;
    
    %Energy per bit to noise density
    EbN0=4.5;
    R_real_i=10^(-0.1*(EbN0-P_rx_i+N0));

    R_real=[R_real; R_real_i];
    r=[r; r_i_norm];
end

XDates = linspace(datetime(2011,3,18), datetime(2015,4,30), length(date));

figure()
yyaxis left
plot(XDates, r/149597870700, '--', 'Color', [0 0.4470 0.7410]  ,'LineWidth', 1)
xlabel('Date')
ylabel('Distance [AU]') 

yyaxis right
plot(XDates, R_real, 'Color', [0.6350 0.0780 0.1840],  'LineWidth', 1)
xlabel('Date')
ylabel('R_{real} [bps]')


legend('Earth-Mercury Distance','Data Rate', 'Location','northwest','Box','on')
grid on
title('Earth-Mercury distance and Data Rate for LGA')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   RAGGIO - HIGH GAIN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%MODULATION AND ENCODING
alpha_RS=255/223;
alpha_conv=6;
alpha_enc= alpha_conv*alpha_RS;
alpha_mod=1; 

f=8.4e9; %Hz
P_in=42; %76.6; 

mu_amp=0.32; %from graph
BER=1e-5; %literature

P_tx=mu_amp*P_in;
P_tx=10*log10(P_tx);

%GAIN
c=299792458;
lambda=c/f;
G_ph=10*log10(8) + 17;


L_cable= -1; %dB

%GROUND
Drx=34;
mu_rx=0.55;
Grx=10*log10(pi^2*Drx^2*mu_rx/lambda^2);
theta_rx=65.3*lambda/Drx;

    %point losses
L_point=-0.1; %dB DSN datasheet
etarx=sqrt(L_point/(-12))*theta_rx;
    
    %Atmospheric loss 
    L_atm=-0.80; %dB from graph 
    
    % Effective Isotropic Radiated Power
    EIRP=P_tx+G_ph+L_cable;

    Kb=1.38e-23; %boltzman
    Ts = 21; %Kelvin
    N0 = 10*log10(Kb*Ts);

mu_S=astroConstants(4);

% Initial date
date_i=date2mjd2000([2011 03 18 0 0 0]);

% Final date
date_f=date2mjd2000([2015 04 30 0 0 0]);

date=date_i:1:date_f;

R_real=[];
r=[];
SNR_B=[];

for i=1:length(date)
    

    % Mercury position
    [kep_m,ksun_m] = uplanet(date(i), 1);
    [r_m,v_m] = kep2car(kep_m(1), kep_m(2), kep_m(3), kep_m(4), kep_m(5), kep_m(6), mu_S);

    % Earth position
    [kep_e,ksun_e] = uplanet(date(i), 3);
    [r_e,v_e] = kep2car(kep_e(1), kep_e(2), kep_e(3), kep_e(4), kep_e(5), kep_e(6), mu_S);

    % Earth - Mercury distance
    r_i=r_m-r_e;
    r_i_norm=norm(r_i)*1e3;

    % Space loss
    L_space=20*log10(lambda/(4*pi*r_i_norm));

    % Received power
    P_rx_i = EIRP + Grx + L_space + L_point + L_atm;
    
    %Energy per bit to noise density
    EbN0=4.5;
    R_real_i=10^(-0.1*(EbN0-P_rx_i+N0));

    R_real=[R_real; R_real_i];
    r=[r; r_i_norm];
end

XDates = linspace(datetime(2011,3,18), datetime(2015,4,30), length(date));

figure()
yyaxis left
plot(XDates, r/149597870700, '--', 'Color', [0 0.4470 0.7410],'LineWidth', 1)
xlabel('Date')
ylabel('Distance [AU]')

yyaxis right
plot(XDates, R_real, 'Color', [0.6350 0.0780 0.1840],  'LineWidth', 1)
xlabel('Date')
ylabel('R_{real} [bps]')


legend( 'Earth-Mercury Distance','Data Rate','Location','northwest','Box','on')
grid on
title('Earth-Mercury distance and Data Rate for PAA')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   RAGGIO - Medium GAIN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%MODULATION AND ENCODING
alpha_RS=255/223;
alpha_conv=6;
alpha_enc= alpha_conv*alpha_RS;
alpha_mod=1; 


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

   %point losses
L_point=-0.1; %dB DSN datasheet
etarx=sqrt(L_point/(-12))*theta_rx;
    
    %Atmospheric loss 
    L_atm=-0.80; %dB from graph
    
    % Effective Isotropic Radiated Power
    EIRP=P_tx+G_ph+L_cable;

    Kb=1.38e-23; %boltzman
    Ts = 21; %Kelvin
    N0 = 10*log10(Kb*Ts);

mu_S=astroConstants(4);

% Initial date
date_i=date2mjd2000([2011 03 18 0 0 0]);

% Final date
date_f=date2mjd2000([2015 04 30 0 0 0]);

date=date_i:1:date_f;

R_real=[];
r=[];
SNR=[];

for i=1:length(date)
    

    % Mercury position
    [kep_m,ksun_m] = uplanet(date(i), 1);
    [r_m,v_m] = kep2car(kep_m(1), kep_m(2), kep_m(3), kep_m(4), kep_m(5), kep_m(6), mu_S);

    % Earth position
    [kep_e,ksun_e] = uplanet(date(i), 3);
    [r_e,v_e] = kep2car(kep_e(1), kep_e(2), kep_e(3), kep_e(4), kep_e(5), kep_e(6), mu_S);

    % Earth - Mercury distance
    r_i=r_m-r_e;
    r_i_norm=norm(r_i)*1e3;

    % Space loss
    L_space=20*log10(lambda/(4*pi*r_i_norm));

    % Received power
    P_rx_i = EIRP + Grx + L_space + L_point + L_atm;
    
    %Energy per bit to noise density
    EbN0=4.5;
    R_real_i=10^(-0.1*(EbN0-P_rx_i+N0));

    R_real=[R_real; R_real_i];
    r=[r; r_i_norm];
end

XDates = linspace(datetime(2011,3,18), datetime(2015,4,30), length(date));

figure()
yyaxis left
plot(XDates, r/149597870700, '--', 'Color', [0 0.4470 0.7410],'LineWidth', 1)
xlabel('Date')
ylabel('Distance [AU]')

yyaxis right
plot(XDates, R_real, 'Color', [0.6350 0.0780 0.1840],  'LineWidth', 1)
xlabel('Date')
ylabel('R_{real} [bps]')


legend( 'Earth-Mercury Distance','Data Rate', 'Location','northwest','Box','on')
grid on
title('Earth-Mercury distance and Data Rate for MGA')
