%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     PROCESSING PARAMETERS                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
warning off

%% PARAMETERS

WorkDir = uigetdir('C:\Users\noe.monnier\Documents\Turbulent analysis\'); % Working directory = directory where images are stored
LaserDir = uigetdir('C:\Users\noe.monnier\Documents\Turbulent analysis\'); % Directory with lasersheet images

% Camera settings

F = 1000; % Camera framerate [fps]
Magn = 0.1716; % Magnification of the optical setup [mm/pxl]
sizex = 512; % Image width [pxl]
sizey = 464; % Image height [pxl]
HIT_center = [256 216]; % HIT zone center [pxl]

% Gas properties

% 0.95 NH3 0.05 H2 ER 1 423K 1bar
rho_b = 0.139285445/1e6; % Density of burnt gas [g/mm^3] 
rho_u = 0.741066152/1e6; % Density of unbrunt gas [g/mm^3]
sL_0 = 0.154343863; % Laminar flame speed

fig_title = "0.95 NH$_3$/0.05 H$_2$/Air $\varphi$=1.0 P=1 bar T=423 K u'=0.49 m/s";

save([WorkDir '\Parameters'], 'F','Magn','sizex','sizey','HIT_center','rho_b','rho_u','sL_0','fig_title','LaserDir') % Save the processing parameters

disp('#  PARAMETERS SAVED  #')