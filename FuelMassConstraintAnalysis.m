%% Fuel vs W/S for given range (includes booster and warhead fractions)
clear
clc
close all
set(groot, 'defaultTextInterpreter','latex')

%% GIVEN (from user)
W_TO = 2000; % lb (takeoff weight)

booster_frac = 0.17; % fraction of takeoff weight (booster mass fraction)
W_booster = W_TO * booster_frac; % lb
W_cruise = W_TO - W_booster; % lb (weight entering cruise phase)

warhead_frac = 0.3; % fraction of cruise weight that is warhead
W_warhead = W_cruise * warhead_frac; % lb

WTO_S_range = 200:50:1000; % lb/ft^2

rho = 8.91e-4;          % density at cruise altitude (slugs/ft^3)
V = 3.5 * 678.1 * 1.46667; % cruise velocity ft/s (Mach * speed of sound in ft/s)
CD0 = 0.12;             % zero-lift drag coefficient
e = 0.7;                % span efficiency
AR = 7^2 / (S_X51);     % example aspect ratio, or calculated from geometry
c_T = 0.54 / 3600; % TSFC (1/s)
R = 500*6076; %Range (ft)


% Step 1: wing area
S = W_warhead / WTO_S_range; % ft^2

% Step 2: dynamic pressure
q = 0.5 * rho * V^2;    % lb/ft^2

% Step 3: lift coefficient at cruise
CL = WTO_S_range ./ q;

% Step 4: induced drag factor
k = 1 / (pi * e * AR);

% Step 5: total drag coefficient
CD = CD0 + k * CL^2;

LD = CL./CD;

% Step 6: fuel mass via Breguet
weight_ratio = exp(- (c_T * R) / (V * LD));
W_final = W_cruise * weight_ratio;
fuel_mass = W_cruise - W_final; % lb



