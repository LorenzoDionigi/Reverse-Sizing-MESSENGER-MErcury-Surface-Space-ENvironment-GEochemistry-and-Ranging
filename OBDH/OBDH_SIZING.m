%OBDH REVERSE SIZING
%
% FIRST OPERATING MISSION DAY AROUND MERCURY 18/03/2011
% LAST OPERATING MISSION DAY AROUND MERCURY 30/04/2015
% 32-bit logic
clc
clear all
close all 

% MEMORY SIZING
% GENERAL CODE
Kbytes_NP=20; % Nucleus plus #kbytes
Kbytes_Clf = 200; % C library functions #kbytes
Kbytes_AP = 825; % Applications program #kbytes
Kbytes_S =200; %Stack #kbytes
Kbytes_SC = 25; %Stack catalog #kbytes
Kbytes_Ars = 20; %Autonomy rule storage #kbytes
Kbytes_Sc = 80; %Stored Commands #kbytes

Total_Memory_Code = (Kbytes_NP+Kbytes_Clf+Kbytes_AP+Kbytes_S+Kbytes_SC+Kbytes_Ars+Kbytes_Sc); %kB
%total is 1370Kbytes (fonte MESSENGER)
Total_Memory_Code = Total_Memory_Code/125; %Mb

%INSTRUMENTS
kB_MDIS = 65; % Memory usage kB
kB_GRNS = 30; % Memory usage kB
kB_MAG = 30; % Memory usage kB
kB_MLA = 20; % Memory usage kB
kB_ASCS = 20; % Memory usage kB
kB_EPPS = 50; % Memory usage kB
kB_XRS = 30; % Memory usage kB
kB_Common_sw = 30; % Memory usage kB

Total_Memory_Instr_Code = (kB_MDIS+kB_GRNS+kB_MAG+kB_MLA+kB_ASCS+kB_EPPS+kB_XRS+kB_Common_sw); %kB
%total is 275 (54% of margin? COSA è STA PERCENTUALE?)
Total_Memory_Instr_Code = Total_Memory_Instr_Code/125; %Mb

num_day = 365.25*4+33;
% mostrare qunato serve per storare i data per 1 mese 6 mesi
% 1 anno e tutta la missione 4 anni +33 giorni e far vedere
% che in caso di perdita di comunicazione con la terra si possono storare
% dati per circa 6 mesi.

%Daily data volume
DV_MDIS=15; %Mb/day
DV_GRNS=6.9; %Mb/day
DV_MAG=0.5; %Mb/day
DV_MLA=2.7; %Mb/day
DV_ASCS=5.4; %Mb/day
DV_EPPS=6.8; %Mb/day
DV_XRS=3.4; %Mb/day
 
Total_Data_volume = (DV_MDIS+DV_GRNS+DV_MAG+DV_MLA+DV_ASCS+DV_EPPS+DV_XRS)*num_day; %Mb

Minimum_Memory_Size=Total_Memory_Instr_Code+Total_Memory_Code+Total_Data_volume
%is it lower than 8Gbits?
Minimum_Memory_Size= Minimum_Memory_Size/1000 %Gb

%OBC features as frequency and throughputs by similiarit
%frequency is 25 or 6.25 Mhz in base alla mode.
%la frequenza del FPP è di 10Mhz
%per i dati OBC commentare che le percentuali richieste dal 
%processore riportate in tabella MP CPU Usage

