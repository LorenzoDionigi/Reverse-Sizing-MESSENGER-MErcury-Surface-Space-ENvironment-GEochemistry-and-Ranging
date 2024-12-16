%% ATTITUDE DETERMINATION and CONTROL SYSTEM (ADCS) of MESSENGER (3-axis stabilized s/c)
clc
clear all
close all

%% DATA
mu_M = astroConstants(11); % Mercury planetary constant [km^3/s^2]
mu_S = astroConstants(4);  % Sun planetary constant      [km^3/s^2]
m = 1107.9;                % Messenger wet mass [Kg] {NASA: https://solarsystem.nasa.gov/missions/messenger/in-depth/}
point_acc = 0.1;           % Messenger pointing accuracy FOR EACH CONTROL MODE [deg] {MESSENGER del 99, page 108 pdf}
% theta                    % Maximum deviation of z-axis from the local vertical [rad]
I_max = 533.5;             % Maximum moment of inertia (around z-axis) [Kg*m^2] {Telemetry Recovery and Uplink Commanding of a Spacecraft Prior to Three-Axis Attitude Stabilization}
I_med=449.6;               % Medium moment of inertia (around x-axis) [Kg*m^2] {Telemetry Recovery and Uplink Commanding of a Spacecraft Prior to Three-Axis Attitude Stabilization
I_min = 432.5;             % Minimum moment of inertia (around y-axis) [Kg*m^2] {Telemetry Recovery and Uplink Commanding of a Spacecraft Prior to Three-Axis Attitude Stabilization}
a_mess = 10179.2497;            % Semi-major axis around Mercury [Km]
e_mess = 0.7399765;             % Eccentricity around Mercury [-]
                           % (a and e are retreived, as a first approximation, from ephemerides
                           % of Messenger orbit after MOI (performed in reality at 00:45 on March
                           % 18th 2011) at 01:15:01 on March 18th 2011. Options: target body
                           % Messenger (spacecraft), coordinate center Mercury (body center) [500@199])
I = 0;                     % Incidence angle in worst case [deg]
% d = c_sp - c_g           % Distance center of solar pressure - center of gravity [m]
th_dot_max_nom = 0.8;      % Maximum slew rate [deg/s] {MESSENGER del 99, page 108 pdf}
h_max_RW = 7.5;            % Maximum angular momentum (storage) of the RW [Nms] {Vaughan, The Messenger Spacecraft Guidance and Control System}
T_max_RW = 0.075;          % Maximum torque of the RW [Nm] {Vaughan, The Messenger Spacecraft Guidance and Control System}
I_s_44 = 227;              % 4.4N-thruster specific impulse [s]
I_s_22 = 234;              % 22N-thruster specific impulse [s]
F_44 = 4.4;                % Thrust of 4.4N-thruster (max) [N]
F_22 = 22;                 % Thrust of 22N-thruster (max) [N]
n_44 = 12;                 % Number of 4.4N-thrusters [number of thrusters] {Homework2}
n_22 = 4;                  % Number of 22N-thrusters [number of thrusters] {Homework2}
% L                        % Arm of EACH thruster [m]
% tb                       % Burning time of EACH TYPE of thruster [s]

% omega_x0 = 0.01; % Tip-off rotation before de-tumbling [rad/s]
% omega_y0 = 0.01; % Tip-off rotation before de-tumbling [rad/s]
% omega_z0 = 0.03; % Tip-off rotation before de-tumbling [rad/s]

% %% Statistical formulas
% % For preliminary sizing we exploit the mass value to understand wich are
% % the best sensors and actuators
% V = 0.01*m;            % Messenger volume [m^3]
% L = 0.25*(m^(1/3));    % Messenger length [m]
% A = L^2;               % Messenger area [m^2]
% I_sc = 0.01*(m^(5/3)); % Messenger moment of inertia [Kgm^2]

%% External disturbances 

% ------------------------ Gravity gradient -------------------------------
%b = a*sqrt(1-(e^2));                                    % Semi-minor axis around Mercury [Km]
%R = sqrt(a*b);                                          % Radius of a circular orbit with the same area of the
                                                        % real one (the elliptical one) [Km]
%T_gg = ((3*mu_M)/(2*(R^3)))*(I_max-I_min)*sin(2*theta); % Gravity Gradient disturbance torque [Nm]

% ------------------------------ SRP --------------------------------------
%Fs = P*c;                           % Solar constant around Mercury [W/m^2] {Lab07-SRP and magnetic disturbance torque, SAD course}
%T_srp = (Fs/c)*A_sh*(1+q)*cos(I)*d; % SRP disturbance torque (on the sunshade?) [Nm]

%% Internal disturbances
% c_g_unc = 0.03;     % Uncertainty in c_g (worst case) [m]
% thruster_mis = 0.5; % Thruster misalignment (worst case) [deg]
% % thruster_mism = +- 5%
%% SLEW
%Slew around z-axis for SKI:
th_max2 = 20;                               % Slew angle around z-axis for SKI slew [deg]
t_min2 = th_max2/th_dot_max_nom;            % SKI slew minimum time [s]
T_slew_min2 = (4*th_max2*I_max/(t_min2^2)); % Torque around z-axis for SKI slew [Nm]
b=sqrt(3)/4;
B=[-b -b b; b -b b; b b b; -b b b];         % Pseudo-inverse matrix
T_rw=B*[0; 0; T_slew_min2];                 % Torque vector (each position is a reaction wheel) [Nm]
if abs(T_rw(:)) > T_max_RW
    fprintf('Slew maneuver not possible. An angular velocity decrease is applied.');
    a=1/sqrt(3);
    A=[-a a a -a; -a -a a a; a a a a];
    T_slew_max2 = A*T_max_RW*ones(4,1); % Maximum torque around z-axis for SKI slew [Nm]
    t_min2 = sqrt(4*th_max2*I_max/T_slew_max2(3)); % Slew minimum time with maximum RW torque [s]
    th_dot_max2 = th_max2/t_min2; % Maximum slew rate with maximum RW torque [deg/s]
    t_min2_minutes=t_min2/60; %Minimum allowed slew time in minutes [min]
    if t_min2>30*60
        error('The maneuver takes more than 30 min')
    else
        disp('The maneuver is possible in less than 30 min')
    end
else
    fprintf("Slew maneuver is possible at the highest angular velocity")
end
%%
% Slew around x-axis for SKI:
th_max3 = 24;                               % Slew angle around x-axis for SKI slew [deg]
t_min3 = th_max3/th_dot_max_nom;            % SKI slew minimum time [s]
T_slew_min3 = (4*th_max3*I_med/(t_min3^2)) % Torque around x-axis for SKI slew [Nm]
b=sqrt(3)/4;
B=[-b -b b; b -b b; b b b; -b b b];         % Pseudo-inverse matrix
T_rw=B*[T_slew_min3; 0; 0];                 % Torque vector (each position is a reaction wheel) [Nm]
if abs(T_rw(:)) > T_max_RW
    fprintf('Slew maneuver not possible. An angular velocity decrease is applied.');
    a=1/sqrt(3);
    A=[-a a a -a; -a -a a a; a a a a];
    T_slew_max3 = A*T_max_RW*[-1;1;1;-1]; % Maximum torque around x-axis for SKI slew [Nm]
    t_min3 = sqrt(4*th_max3*I_med/T_slew_max3(1)); % Slew minimum time with maximum RW torque [s]
    th_dot_max3 = th_max3/t_min3; % Maximum slew rate with maximum RW torque [deg/s]
    t_min3_minutes=t_min3/60; %Minimum allowed slew time in minutes [min]
    if t_min3>30*60
        error('The maneuver takes more than 30 min')
    else
        disp('The maneuver is possible in less than 30 min')
    end
else
    fprintf("Slew maneuver is possible at the highest angular velocity")
end

%% WORST CASE SLEW SCENARIO FOR Y AXIS
th_max4 = 180;                              % Slew angle around y-axis for worstcase slew [deg]
t_min4 = th_max4/th_dot_max_nom;            % Slew minimum time [s]
T_slew_min4 = (4*th_max4*I_min/(t_min4^2)); % Torque around y-axis for worstcase slew [Nm]
b=sqrt(3)/4;
B=[-b -b b; b -b b; b b b; -b b b];         % Pseudo-inverse matrix
T_rw=B*[0; T_slew_min4; 0];                 % Torque vector (each position is a reaction wheel) [Nm]
if abs(T_rw(:)) > T_max_RW
    disp('Slew maneuver not possible. An angular velocity decrease is applied.');
    a=1/sqrt(3);
    A=[-a a a -a; -a -a a a; a a a a];
    T_slew_max4 = A*T_max_RW*[-1; -1; 1;1]; % Maximum torque around y-axis for worstcase slew [Nm]
    t_min4 = sqrt(4*th_max4*I_min/T_slew_max4(2)); % Slew minimum time with maximum RW torque [s]
    th_dot_max4 = th_max4/t_min4; % Maximum slew rate with maximum RW torque [deg/s]
    t_min4_minutes=t_min4/60; %Minimum allowed slew time in minutes [min]
    if t_min4>30*60
        error('The maneuver takes more than 30 min')
    else
        fprintf('The maneuver is possible, in less than 30 min')
    end
else
    disp("Slew maneuver is possible at the highest angular velocity")
end

%% DISTURBANCES: SRP
% DATA:
Fs = 1358; % Solar constant at 1 AU [W/m^2] {Lab Bernelli}
AU = astroConstants(2); % Astronomical unit [Km]
cs_sh = 0.15;              % Sunshade specular reflectivity coefficient [-] {Vaughan, MOMENTUM MANAGEMENT FOR THE MESSENGER MISSION, page 14 pdf}
cd_sh = 0.35;              % Sunshade diffuse reflectivity coefficient [-] {Vaughan, MOMENTUM MANAGEMENT FOR THE MESSENGER MISSION, page 14 pdf}
cs_sp = 0.52;              % Solar panel specular reflectivity coefficient [-] {Vaughan, MOMENTUM MANAGEMENT FOR THE MESSENGER MISSION, page 14 pdf}
cd_sp = 0.07;              % Solar panel diffuse reflectivity coefficient [-] {Vaughan, MOMENTUM MANAGEMENT FOR THE MESSENGER MISSION, page 14 pdf}
%P = 5.0145e-5;             % Solar pressure (at 0.3 AU distance from the Sun) [N/m^2] {Vaughan, MOMENTUM MANAGEMENT FOR THE MESSENGER MISSION, page 14 pdf}
c = 299792458;             % Speed of light [m/s]
A_sh = 5.24553;            % Sunshade surface toward the Sun [m^2] {Berman, Domingue, et al, Orbital Operations Planning and Scheduling for the MESSENGER Mission, APL}
                           % (This is the semi-cylindrical area of the
                           % sunshade but, to simplify we assume that it is
                           % a flat plate)
A_sp = 5;                  % Solar panels surface toward the Sun (2.5 each) [m^2] {Vaughan, MOMENTUM MANAGEMENT FOR THE MESSENGER MISSION, page 14 pdf}

%% Keplerian parameters from MOI (2 days)
% % Remember to import data (C589:N20068) from Ephemerides_2days_Mercury
% a_v = Ephemerides2daysMercury(:,10); % [Km]
% e_v = Ephemerides2daysMercury(:,1); % [-]
% i_v = Ephemerides2daysMercury(:,3); % [deg]
% OM_v = Ephemerides2daysMercury(:,4); % [deg]
% om_v = Ephemerides2daysMercury(:,5); % [deg]
% th_v = Ephemerides2daysMercury(:,9); % [deg]
%% Keplerian parameters from MOI (88 days)
% Remember to import data (C80:L20068) from Ephemerides_88days_Mercury
a_v = Ephemerides88daysMercury(:,10); % [Km]
e_v = Ephemerides88daysMercury(:,1); % [-]
i_v = Ephemerides88daysMercury(:,3); % [deg]
OM_v = Ephemerides88daysMercury(:,4); % [deg]
om_v = Ephemerides88daysMercury(:,5); % [deg]
th_v = Ephemerides88daysMercury(:,9); % [deg]
%% Keplerian parameters from MOI of Mercury (2 days)
% % Remember to import data (C570:L20049) from Ephemerides_2days_MercuryBody
% aM_v = Ephemerides2daysMercuryBody(:,10); % [Km]
% eM_v = Ephemerides2daysMercuryBody(:,1); % [-]
% iM_v = Ephemerides2daysMercuryBody(:,3); % [deg]
% OMM_v = Ephemerides2daysMercuryBody(:,4); % [deg]
% omM_v = Ephemerides2daysMercuryBody(:,5); % [deg]
% thM_v = Ephemerides2daysMercuryBody(:,9); % [deg]
%% Distance from the Sun and Solar pressure
R_mat = [];
D_s_vect = [];
P_vect = [];
for cont = 1:length(a_v)
    [r_ins,v_ins] = kep2car(a_v(cont,1), e_v(cont,1), deg2rad(i_v(cont,1)),...
        deg2rad(OM_v(cont,1)), deg2rad(om_v(cont,1)), deg2rad(th_v(cont,1)), mu_S);
    R_mat = [R_mat; r_ins'];
    D_s = norm(r_ins);
    D_s = D_s/AU; % Distance from the Sun [AU]
    P = Fs/(c*(D_s^2)); % Solar pressure at the distance D_s [Pa?]
    D_s_vect = [D_s_vect D_s]; % Distances from the Sun [AU]
    P_vect = [P_vect P]; % Solar pressure along 2 days (4 orbits around Mercury) [Pa?]
end
%%
figure(1)
plot(1:length(D_s_vect),D_s_vect)
title('Distance from the Sun after MOI along 2 days')
figure(2)
plot3(R_mat(:,1),R_mat(:,2),R_mat(:,3))
%% Integration of Fs (for EPS homework)
Fs_vect = [];
for cont = 1:length(a_v)
    Fs_ = P_vect(cont)*c; % Solar pressure at the distance D_s [Pa?] P_vect(cont)*(c*(D_s_vect(cont)^2))
    Fs_vect = [Fs_vect Fs_];
end
Fs_integral = trapz(D_s_vect,Fs_vect);
Fs_int = trapz(Fs_vect)/(max(D_s_vect)-min(D_s_vect));
plot(D_s_vect,Fs_vect)
%% Function for theta solar panels

% % If 2 days are considered:
% D_s_chart = [0.30 0.323 0.35 0.367 0.40 0.413 0.423 0.44 0.447 0.455 0.458 0.46]; % Solar distance from chart [AU] 0.45
% th_sp_chart = [64.5 60 54 50 40 35 30 20 15 5 2.5 0]; % Solar array angle constraint from chart [deg] 10
% span = 0.0001;
% D_s_plot1 = 0.30:span:(min(D_s_vect)-span);
% D_s_plot2 = (max(D_s_vect)+span):span:0.46; 
% D_s_plot = [D_s_plot1 D_s_vect D_s_plot2]; % Vector of solar distances to evaluate the solar array angle constraint [AU]

% If 88 days are considered:
D_s_chart = [0.30 0.323 0.35 0.367 0.40 0.413 0.423 0.44 0.447 0.455 0.458 0.46]; % Solar distance from chart [AU] 0.45
th_sp_chart = [64.5 60 54 50 40 35 30 20 15 5 2.5 0]; % Solar array angle constraint from chart [deg] 10
span = 0.0001;
D_s_plot1 = 0.30:span:(min(D_s_vect)-span);
D_s_vect_lim = [];
i_lim1 = 1;
D_s_vect = D_s_vect*AU; % [Km]
[D_s_max,i_max] = max(D_s_vect);
D_s_vect = D_s_vect(1:i_max); % [Km]
P_vect = P_vect(1:i_max);
for i_plot=1:length(D_s_vect) % Only distances under 0.46 are taken to apply a fitting process (above 0.46 AU, th_ins = 0)
    if D_s_vect(i_plot) <= (0.46*AU)
        D_s_vect_lim(i_lim1) = D_s_vect(i_plot)/AU; % [AU]
        i_lim1 = i_lim1 + 1;
    end
end
D_s_plot = [D_s_plot1 D_s_vect_lim]; % Vector of Solar distances for the chart of next section [AU]
D_s_vect = D_s_vect/AU; % [AU]

figure(9)
plot(1:length(D_s_vect),D_s_vect)
title('Distance from the Sun after MOI along 44 days')
%%
figure(10)
plot(D_s_vect,P_vect)
title('Solar pressure after MOI along 44 days')
xlabel('Solar distance [AU]'); ylabel('Solar pressure [N/m^2]')
%% Interpolation for Solar panel angle vs Solar distance function
% Polyfit
th_sp_pol1 = polyfit(D_s_chart,th_sp_chart,2);
th_sp1 = polyval(th_sp_pol1,D_s_plot);
% Spline
th_sp2 = spline(D_s_chart,th_sp_chart,D_s_plot);
% Interp1
th_sp3 = interp1(D_s_chart,th_sp_chart,D_s_plot);
%%
figure(3)
plot(D_s_plot,th_sp1,'-b',D_s_plot,th_sp2,'-r',D_s_plot,th_sp3,'-k','LineWidth',1.5); %---> interp1 better approx
title('Solar array angle constraint')
legend('Polyfit','Spline','Interp1')
xlabel('Solar distance [AU]'); ylabel('Solar array rotation angle [deg]');
grid on
figure(4)
plot(D_s_plot,th_sp2,'-r','LineWidth',1.5); 
title('Solar array angle constraint')
xlabel('Solar distance [AU]'); ylabel('Solar array rotation angle [deg]');
grid on
figure(5)
plot(D_s_plot,th_sp3,'-k','LineWidth',1.5);
title('Solar array angle constraint')
xlabel('Solar distance [AU]'); ylabel('Solar array rotation angle [deg]');
grid on

%% SRP torque computation

TSRP_x12_vect = [];
TSRP_xsh_vect = [];
TSRP_x_vect = [];
F_sp_y12_vect = [];
F_sp_z12_vect = [];
F_sh_y_vect = [];
A_sp_cross_vect = [];
for cont = 1:length(D_s_vect)
    D_s = D_s_vect(cont);
    P = P_vect(cont);
    for cont3 = 1:length(D_s_plot)
        if D_s == D_s_plot(cont3)
            th_ins = th_sp3(cont3);
        end
        if D_s > 0.46
            th_ins = 0;
        end
    end
    A_sp_cross = (A_sp/2)*cos(deg2rad(th_ins)); % Cross-sectional area of one solar panel [m^2]
    A_sp_cross_vect = [A_sp_cross_vect A_sp_cross];

    % Assumption: center of solar pressure in the geometric center of the
    % solar panel.
    c_sp1_x = 2.28; % Center location along x-axis of the solar panel [m]

    c_sp12_y = -0.008; % y component of spacecraft center of mass [m]
    c_sp12_z = -0.055; % z component of spacecraft center of mass [m]

    %b_sp1_x = c_sp1_x;    % Arm of the SRP (along x-axis) acting on the solar panel [m]
    b_sp12_y = -c_sp12_y; % Arm along y of the forces along z [m]
    b_sp12_z = c_sp12_z;  % Arm along z of the forces along y [m]

    offset_spc_cm = 0.05; % Possible offset between the shade center and the spacecraft center of mass, along z [m]

    % Forces:
    F_sp_y12 = -P*A_sp_cross*(-(1-cs_sp) - 2*( cs_sp*cos(deg2rad(th_ins)) + (1/3)*cd_sp )*cos(deg2rad(th_ins)));
    F_sp_z12 = -P*A_sp_cross*(-2*( cs_sp*cos(deg2rad(th_ins)) + (1/3)*cd_sp )*sin(deg2rad(th_ins)));
    F_sh_y = P*A_sh*(1 + cs_sh + (2/3)*cd_sh);
    F_sp_y12_vect = [F_sp_y12_vect F_sp_y12];
    F_sp_z12_vect = [F_sp_z12_vect F_sp_z12];
    F_sh_y_vect = [F_sh_y_vect F_sh_y];
    % Torques:
    % - due to F_sp_y12
          %MSRP_z1_y = F_sp_y12*b_sp1_x; % si annulla con l'altro pannello
          TSRP_x1_y = F_sp_y12*b_sp12_z; % for one sp
          TSRP_x12_y = TSRP_x1_y*2; % for 2 sp
    % - due to F_sp_z12
          %MSRP_y1_z = F_sp_z12*b_sp1_x; % si annulla con l'altro pannello
          TSRP_x1_z = F_sp_z12*b_sp12_y; % for one sp
          TSRP_x12_z = TSRP_x1_z*2; % for 2 sp
    TSRP_x12 = TSRP_x12_y + TSRP_x12_z; % total around x for the solar panels
    % - due to F_sh_y
          TSRP_xsh = F_sh_y*offset_spc_cm;
          if TSRP_xsh*TSRP_x12 < 0
              TSRP_xsh = -TSRP_xsh;
          end
    TSRP_x = TSRP_xsh + TSRP_x12;
    TSRP_x12 = norm(TSRP_x12);
    TSRP_xsh = norm(TSRP_xsh);
    TSRP_x12_vect = [TSRP_x12_vect TSRP_x12];
    TSRP_xsh_vect = [TSRP_xsh_vect TSRP_xsh];
    TSRP_x = norm(TSRP_x);
    TSRP_x_vect = [TSRP_x_vect TSRP_x];
end

%% SRP torque plot as a function of distance from the Sun
% (During the first two days after MOI)
figure(6)
plot(D_s_vect,TSRP_x_vect,'-b','LineWidth',2)
title('SRP torque plot vs Distance from the Sun')
xlabel('Distance from the Sun [AU]'); ylabel('SRP Torque around x-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(TSRP_x_vect) max(TSRP_x_vect)])

figure(7)
plot(D_s_vect,TSRP_xsh_vect,'-r','LineWidth',2)
title('Sunshade SRP torque')
xlabel('Distance from the Sun [AU]'); ylabel('SRP Torque around x-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(TSRP_xsh_vect) max(TSRP_xsh_vect)])

figure(8)
plot(D_s_vect,TSRP_x12_vect,'-g','LineWidth',2)
title('Solar panels SRP torque')
xlabel('Distance from the Sun [AU]'); ylabel('SRP Torque around x-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(TSRP_x12_vect) max(TSRP_x12_vect)+(2e-7)])

figure(11)
plot(D_s_vect,F_sp_y12_vect,'-k','LineWidth',2)
title('Solar panel SRP force along y-axis')
xlabel('Distance from the Sun [AU]'); ylabel('SRP force along y-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(F_sp_y12_vect) max(F_sp_y12_vect)+(2e-7)])

figure(12)
plot(D_s_vect,F_sp_z12_vect,'-k','LineWidth',2)
title('Solar panel SRP force along z-axis')
xlabel('Distance from the Sun [AU]'); ylabel('SRP force along z-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(F_sp_z12_vect) max(F_sp_z12_vect)])

figure(13)
plot(D_s_vect,F_sh_y_vect,'-k','LineWidth',2)
title('Sunshade SRP force along y-axis')
xlabel('Distance from the Sun [AU]'); ylabel('SRP force along y-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(F_sh_y_vect) max(F_sh_y_vect)])
%%
figure(14)
plot(D_s_vect,A_sp_cross_vect,'-c','LineWidth',1.5)
title('Cross-sectional area of solar panels')
xlabel('Distance from the Sun [AU]'); ylabel('A_s_p cross-sectional [m^2]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(A_sp_cross_vect) max(A_sp_cross_vect)])

%% tiledlayout
figure(15)
tiledlayout(3,1)
nexttile
plot(D_s_vect,TSRP_x_vect,'-b','LineWidth',2)
title('SRP torque plot vs Distance from the Sun')
xlabel('Distance from the Sun [AU]'); ylabel('SRP Torque around x-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(TSRP_x_vect) max(TSRP_x_vect)])
nexttile
plot(D_s_vect,TSRP_xsh_vect,'-r','LineWidth',2)
title('Sunshade SRP torque')
xlabel('Distance from the Sun [AU]'); ylabel('SRP Torque around x-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(TSRP_xsh_vect) max(TSRP_xsh_vect)])
nexttile
plot(D_s_vect,TSRP_x12_vect,'-g','LineWidth',2)
title('Solar panels SRP torque')
xlabel('Distance from the Sun [AU]'); ylabel('SRP Torque around x-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(TSRP_x12_vect) max(TSRP_x12_vect)+(2e-7)])

%% tiledlayout 
figure(16)
tiledlayout(1,3)

nexttile
plot(D_s_vect,F_sp_y12_vect,'-k','LineWidth',2)
title('Solar panel SRP force along y-axis')
xlabel('Distance from the Sun [AU]'); ylabel('SRP force along y-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(F_sp_y12_vect) max(F_sp_y12_vect)+(2e-7)])
nexttile
plot(D_s_vect,F_sp_z12_vect,'-k','LineWidth',2)
title('Solar panel SRP force along z-axis')
xlabel('Distance from the Sun [AU]'); ylabel('SRP force along z-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(F_sp_z12_vect) max(F_sp_z12_vect)])
nexttile
plot(D_s_vect,F_sh_y_vect,'-k','LineWidth',2)
title('Sunshade SRP force along y-axis')
xlabel('Distance from the Sun [AU]'); ylabel('SRP force along y-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(F_sh_y_vect) max(F_sh_y_vect)])

%%
figure(18)
tiledlayout(1,2)
nexttile
plot(D_s_vect,TSRP_x_vect,'.-b',D_s_vect,TSRP_xsh_vect,'.-r',D_s_vect,TSRP_x12_vect,'-g','LineWidth',1.5)
title('SRP torque vs Distance from the Sun')
xlabel('Distance from the Sun [AU]'); ylabel('SRP Torque around x-axis [Nm]')
legend('Total SRP torque','Sunshade SRP torque','Solar panels SRP torque')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(TSRP_x12_vect) max(TSRP_x_vect)])
nexttile
plot(D_s_vect,F_sp_y12_vect,'.-k',D_s_vect,F_sp_z12_vect,'.-b',D_s_vect,F_sh_y_vect,'.-r','LineWidth',1.5)
title('SRP force')
xlabel('Distance from the Sun [AU]'); ylabel('SRP force [N]')
legend('Solar panels SRP along y','Solar panels SRP along z','Sunshade SRP along y')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(F_sp_z12_vect) max(F_sh_y_vect)+(2e-7)])

%% tiledlayout
figure(17)
tiledlayout(2,3)
nexttile
plot(D_s_vect,TSRP_x_vect,'-b','LineWidth',2)
title('Total SRP torque')
xlabel('Distance from the Sun [AU]'); ylabel('SRP Torque around x-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(TSRP_x_vect) max(TSRP_x_vect)])
nexttile
plot(D_s_vect,TSRP_xsh_vect,'-r','LineWidth',2)
title('Sunshade SRP torque')
xlabel('Distance from the Sun [AU]'); ylabel('SRP Torque around x-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(TSRP_xsh_vect) max(TSRP_xsh_vect)])
nexttile
plot(D_s_vect,TSRP_x12_vect,'-g','LineWidth',2)
title('Solar panels SRP torque')
xlabel('Distance from the Sun [AU]'); ylabel('SRP Torque around x-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(TSRP_x12_vect) max(TSRP_x12_vect)+(2e-7)])
nexttile
plot(D_s_vect,F_sp_y12_vect,'-k','LineWidth',2)
title('Solar panel SRP force along y-axis')
xlabel('Distance from the Sun [AU]'); ylabel('SRP force along y-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(F_sp_y12_vect) max(F_sp_y12_vect)+(2e-7)])
nexttile
plot(D_s_vect,F_sp_z12_vect,'-k','LineWidth',2)
title('Solar panel SRP force along z-axis')
xlabel('Distance from the Sun [AU]'); ylabel('SRP force along z-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(F_sp_z12_vect) max(F_sp_z12_vect)])
nexttile
plot(D_s_vect,F_sh_y_vect,'-k','LineWidth',2)
title('Sunshade SRP force along y-axis')
xlabel('Distance from the Sun [AU]'); ylabel('SRP force along y-axis [Nm]')
grid on
xlim([min(D_s_vect) max(D_s_vect)])
ylim([min(F_sh_y_vect) max(F_sh_y_vect)])
%% ACTUATORS

%% ---------------- REACTION WHEELS (Teldix RSI 7-75/601) ------------------

% External disturbances dump
T_srp = mean(TSRP_x_vect);
T_d = T_srp*2;                         % Overall average disturbance + 100% margin (to be applied with RW)
                                       % at each point around Mercury [Nm]
a_M = 57909070.3576;                   % Mercury semi-major axis [Km] {Ephemerides}
T_orb = 2*pi*sqrt((a_mess^3)/mu_M);    % Orbital period around Mercury [s]
h_d_orb = T_d*T_orb;                   % Average angular momentum accumulated by the RW along 1 orbit [Nms]
n_orb_desat = floor(h_max_RW/h_d_orb); % Number of orbits before RW desaturation [number of orbits]
                                       % (Considering only average disturbances)
%% ---------------------------- THRUSTERS ----------------------------------

% % Generalities-------------------------------------------------------------
% 
% % Torque generate by the thruster system on Messenger
% % T_thruster = n*F*L; OR T_thruster = I_sc*th_2dot; (Messenger momentum inertia and Messenger angular acceleration)
% L_dist = (1.27/2) + 0.052;
% th_dot_max = (n_22*F_22*L)*tb/I_s_22;     % Angular velocity spanned during thrusting phase [rad/s]
% th = 0.5*(n_22*F_22*L)*(tb^2)/I_sc; % Angle spanned during thrusting phase [rad]

%% Disturbances ------------------------------------------------------------
n = 2;
L_dist = (1.27/2) + 0.052;
F_th = T_d/(n*L_dist); % Force to counteract external disturbances [N]

t_dist = 365*24*60*60; % [s] assumption

mp_dist = (t_dist*F_th)/(I_s_22*9.81); % [Kg]
mp_dist_tot = mp_dist*n;

%% Slew before MOI ---------------------------------------------------------

th_dot_MOI = 0.04; % [deg/s] {Guidance and Control Challenges the Messenger}
th_MOI = 30; % Total slew angle [deg]
n = 2; % number of used thrusters
L_22 = (1.27/2) + 0.052;
t_acc_MOI = th_dot_MOI*I_min/(n*F_22*L_22);
th_acc_MOI = 0.5*n*F_22*L_22*(t_acc_MOI^2)/I_min; % Acceleration angle of the slew for MOI, equal to th_brake_MOI [deg]
th_coast_MOI = th_MOI - (2*th_acc_MOI); % Coasting angle [deg]
t_coast_MOI = th_coast_MOI/th_dot_MOI;
tb_MOI = t_coast_MOI + (2*t_acc_MOI);

th_ddot_MOI = th_dot_MOI/(t_acc_MOI*tb_MOI); % Slew acceleration during MOI [deg/s]
F_th_MOI = I_min*th_ddot_MOI/(n*L_22);

tb_22 = sqrt((2*th_MOI*I_min)/(n_22*F_22*L_22));
th_dot_max_22 = (n_22*F_22*L_22)*tb_22/I_min; % Max angular velocity spanned during thrusting phase [rad/s]

% th_dot_max_22 < th_dot_MOI ---> ok, commento su overleaf 

mp_slew = (tb_MOI*F_th_MOI)/(I_s_22*9.81); % [Kg]
mp_slew_tot = mp_slew*n;


%% Desaturation (RW momentum dump) -----------------------------------------

% We could assume n=2 and L=0.5*Messengerlength???
n_desat = 4; % ??
L_44 = (1.27/2) + 0.052; % 1.27 approximation
t_min_desat = h_max_RW/(n_desat*L_44*F_44); % Desaturation minimum time [s]
                                            % Generally, it is around 5 s, then, if
                                            % t_min_desat is too small, take 5 s as tb
% With 4 4.4N-thrusters and L=(1.27/2) + 0.052, we get 0.6203 s ---> take 5s
t_desat = 5; % [s]
F_desat = h_max_RW/(n_desat*L_44*t_desat); % Force to desaturate (one RW?) [N]

% Propellant mass ---------------------------------------------------------
deltat_desat = n_orb_desat*T_orb;
mp_desat = (t_desat*F_desat)/(I_s_44*9.81); % [Kg]
mp_desat_tot = mp_desat*4*(365*24*60*60/deltat_desat); % [Kg] 