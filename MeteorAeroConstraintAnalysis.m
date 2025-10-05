% % constraint analysis for missile design based off of Meteor
clear
clc
close all
set(groot, 'defaultTextInterpreter','latex')

% takeoff, cruise, descent, attack

% because the Meteor is smaller than our requirements, we could
% hypothetically scale it up

% Meteor properties
weight_meteor = 419; % lb
length_meteor = 12; % ft
diameter_meteor = 7/12; % ft (estimation)
range_X51 = 75 * 0.868976; % nm

% very rough approximation of wing-planform area and AR
S_meteor = 4 * 0.75 * 0.65; % ft^2
AR_meteor = 0.75/0.65;
WS_meteor = weight_meteor / S_meteor; % lb/ft^2
TW_meteor_cruise = 0.25; % estimated for ramjet at Mach 4
mach_cruise_meteor = 4;

% our properties based off requirements
W_TO = 2000; % lb
beta_subsonic = 0.7; % estimated for solid rocket
beta_supersonic = 0.85; % estimated for ramjet

W_TO_avg = W_TO * (1 - beta_subsonic/2);
W_TO_end = W_TO * (1 - beta_subsonic/2);

W_cruise_avg = W_TO_end * (1 - beta_supersonic/2);
W_cruise_end = W_TO_end * (1 - beta_supersonic/2);

W_end = W_cruise_end;

mach = 3.5; % based off of what Pat recommended
range = 500; % nm
altitude_cruise = 30000; % ft
altitude_ingress = 500; % low-level ingress
g = 32.2; % ft/s^2
WTO_S_range = 200:50:1000; % lb/ft^2

% assumptions
CL_cruise = 0.15;
CD0_cruise = 0.12;
e = 0.7;
LD_cruise = 8;

% properties at 30,000'
rho_cruise = 8.91e-4; % slugs/ft^3
a_cruise = 678.1 * 1.46667; % ft/s
V_cruise = mach * a_cruise; % ft/s
q_cruise = 0.5 * rho_cruise * V_cruise^2; 

% properties at low-level ingress (sea level), assuming we stay at the same
% speed
V_TO = 1*1116;
rho_ingress = 23.77e-4; % slugs/ft^3
q_ingress = 0.5 * rho_ingress * V_cruise^2;
q_TO = 0.5 * rho_ingress * V_TO^2;

% initiate matrices for TW at different phases of mission
TW_takeoff = zeros(length(WTO_S_range),1);
TW_cruise = zeros(length(WTO_S_range),1);
TW_ingress = zeros(length(WTO_S_range),1);
TW_boost = zeros(length(WTO_S_range),1);


% populate T/W matrices
for i = 1:length(WTO_S_range)
    WS = WTO_S_range(i);
    S_TO = W_TO_avg / WS; % ft^2
    S_cruise = W_cruise_avg / WS; % ft^2

    % Takeoff T/W
    WS_TO_avg = W_TO_avg / S_TO;
    CL_TO = WS_TO_avg / q_TO;
    CD_TO = CD0_cruise + (CL_TO^2 / (pi * AR_meteor * e));
    LD_TO = CL_TO / CD_TO;
    TW_takeoff(i) = 1 / LD_TO;

    % Cruise T/W
    WS_cruise_avg = W_cruise_avg / S_cruise;
    CL_cruise = WS_cruise_avg / q_cruise;
    CD_cruise = CD0_cruise + (CL_cruise^2 / (pi * AR_meteor * e));
    LD_cruise = CL_cruise / CD_cruise;
    TW_cruise(i) = 1 / LD_cruise;

    % Ingress T/W
    WS_ingress = W_end / S_cruise;
    CL_ingress = WS_ingress / q_ingress; % note CL_max = 0.6 constraint can be applied if desired
    CD_ingress = CD0_cruise + (CL_ingress^2 / (pi * AR_meteor * e));
    LD_ingress = CL_ingress / CD_ingress;
    TW_ingress(i) = 1 / LD_ingress;

    % Boost T/W
    accel_target = 80; % ft/s^2 (acceleration requirement during boost)
    TW_boost(i) = accel_target / g;
end

% % analyze T/W & W/S
TW_matrix = [TW_takeoff, TW_cruise, TW_ingress, TW_boost];
TW_max = max(TW_matrix, [], 2);

[min_TW, idx_opt] = min(TW_max);
optimal_WS = WTO_S_range(idx_opt);
optimal_TW = min(TW_max);
%possible = (CL_ingress <= 0.6) & (TW_cruise <= 0.3) & (TW_boost <= 5);
%WTO_S_possible = WTO_S_range(possible);
%TW_max_possible = TW_max(possible);
%if isempty(WTO_S_possible)
%    error('No valid W/S found.')
%end

%[min_TW, i] = min(TW_max_possible);
%optimal_WS = WTO_S_possible(i);

T_needed = optimal_TW*weight_meteor
S_needed = weight_meteor/optimal_WS


% plot T/W vs W/S
plot(WTO_S_range, TW_cruise)
hold on
plot(WTO_S_range, TW_ingress)
plot(WTO_S_range, TW_boost)
plot(WTO_S_range, TW_takeoff)

% plot(optimal_WS, min_TW, '.', 'MarkerSize', 15)
title('Constraint Plot for X-51 Inspired Design')
xlabel('$\frac{W_{TO}}{S} (\frac{lb}{ft^2})$')
ylabel('$\frac{T}{W_{TO}}$')
legend('Cruise','Ingress','Boost','Takeoff')
hold off