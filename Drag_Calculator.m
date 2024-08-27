clc
clear all
close all

%% Given Variables

% Global
rho = 8.754*10^-4;  % Density in slugs/ft^3
T = 400;            % Temperature in degrees Rankine
u = 3.025*10^-7;     % Dynamic Viscosity in (lbs*s)/ft^2
W = 98000;          % Weight in lbs
V_true = 765;       % TAS in ft/s

% Unique Variables

% Fuselage
l_f = 114;          % Fuselage length in feet
d_f = 11;           % Fuselage diameter in feet
Swet_f = 0.8*pi*d_f*l_f;
ratio_f = l_f/d_f;

% Wing
b_wing = 94;        % Wingspan in feet
c_r_wing = 19;           % Centerline root chord in feet
c_r_exp_wing = 17.40;    % Exposed root chord in feet
taper_wing = 0.28; % Wing Taper Ratio
c_t_wing = c_r_wing*taper_wing;
taper_wing_exp = c_t_wing/c_r_exp_wing;      % Exposed Taper Ratio
sweep_wing = 24;      % Quarter Chord Wing Sweep in Degrees
Swing_total = 2361.61; % Total Wing SA in ft^2
Swet_wing = 1943.80;
Sref_wing = 1465.43;   % Planform Wing Area
MAC_wing = 2/3*c_r_exp_wing*(1 + taper_wing_exp-(taper_wing_exp/(1+taper_wing_exp)));   % Mean Aerodynamic Chord in ft
t_MAC_wing = 1.76;       % Max thickness at MAC in ft
tc_MAC_wing = t_MAC_wing/MAC_wing; % Average t/c
AR = b_wing^2/Sref_wing;

% Vertical Tail
S_ref_Vtail = 161;
tc_Vtail = 0.09;
sweep_Vtail = 43.5;
taper_Vtail = 0.8;
c_r_Vtail = 15.5;
MAC_Vtail = 2/3*c_r_Vtail*(1+taper_Vtail-taper_Vtail/(1+taper_Vtail));

% Horizontal Tail
S_ref_Htail = 261;
tc_Htail = 0.09;
sweep_Htail = 31.6;
taper_Htail = 0.35;
c_r_Htail = 11.1;
MAC_Htail = 2/3*c_r_Htail*(1+taper_Htail-taper_Htail/(1+taper_Htail));

% Pylons
S_ref_pylons = 117;
tc_pylons = 0.06;
sweep_pylons = 0;
taper_pylons = 1;
c_pylons = 16.2;
MAC_pylons = 2/3*c_pylons*(1+taper_pylons-taper_pylons/(1+taper_pylons));

% Nacelles
S_wet_n = 455;
ratio_n = 5;
l_n = 16.8;


% Calculations

%% Plot Digitization

% Typical Aircraft Roughness
load("Re.mat")
load("Coeff.mat")

load("Aircraft_Roughness.mat")

Re = linspace(1*10^5,1*10^9,1000);

Aircraft_Roughness = 213.5.*Re.^-0.2841 + 1.082;    % Obtained from Curve Fitter Tool

figure(1)
loglog(Re,Aircraft_Roughness,'-k')
title('Typical Aircraft Surface Roughness')
xlabel('Reynolds Number Based on Length')
ylabel('Skin Friction Coefficient C_f x 10^3')
grid on

% Aero Surface Form Factor
M_0 = 0.8;
tc = linspace(0,0.20,1000);
Z = ((2-M_0^2)*cosd(sweep_wing))/sqrt(1-M_0^2*(cosd(sweep_wing))^2);
K1 = 1 + Z.*tc + 100.*tc.^4;

figure(2)
plot(tc,K1,'-r')
title('Aero Surface Form Factor for Sweep = 24 deg')
xlabel('t/c ratio')
ylabel('Form Factor K')
grid on

% Fineness Ratio

load('K2.mat')
load("L_D.mat")
L_D = linspace(3,11,1000);

K2 = 5.3933*exp(-0.9111*L_D) + 1.3965*exp(-0.0279*L_D);

figure(3)
plot(L_D,K2,'-b')
title('Fineness Ratio vs. Form Factor')
xlabel('Fineness Ratio L/D')
ylabel('Form Factor K')
grid on


%% K and Cf values per Component
v_min = 230;
v_max = 880;
v_step = 0.5;

i_max = (v_max-v_min)/v_step;

% Wing
for i = 1:1:i_max
    v = v_min + (i-1)*v_step;
    Re = (rho*v*MAC_wing)/u;
    wing_roughness = (213.5*Re^-0.2841 + 1.082)/1000;
    Cf_wing(i) = wing_roughness;
end

M_0 = 0.5;
Z = (((2-M_0^2)*cosd(sweep_wing))/sqrt(1-M_0^2*(cosd(sweep_wing))^2))/10;
K1 = 1 + Z*tc_MAC_wing + 100*tc_MAC_wing;

k_wing = K1;

Swet_wing = Swet_wing;

% Fuselage
for i = 1:1:i_max
    v = v_min + (i-1)*v_step;
    Re = (rho*v*l_f)/u;
    fuselage_roughness = (213.5*Re^-0.2841 + 1.082)/1000;
    Cf_f(i) = fuselage_roughness;
end

K2 = 5.3933*exp(-0.9111*ratio_f) + 1.3965*exp(-0.0279*ratio_f);
k_f = K2;

Swet_f = Swet_f;

% Vertical Tail
for i = 1:1:i_max
    v = v_min + (i-1)*v_step;
    Re = (rho*v*MAC_Vtail)/u;
    Vtail_roughness = (213.5*Re^-0.2841 + 1.082)/1000;
    Cf_Vtail(i) = Vtail_roughness;
end

M_0 = 0.5;
Z = (((2-M_0^2)*cosd(sweep_Vtail))/sqrt(1-M_0^2*(cosd(sweep_Vtail))^2))/10;
K1 = 1 + Z*tc_Vtail + 100*tc_Vtail;

k_Vtail = K1;

Swet_Vtail = 1.02*2*S_ref_Vtail;

% Horizontal Tail
for i = 1:1:i_max
    v = v_min + (i-1)*v_step;
    Re = (rho*v*MAC_Htail)/u;
    Htail_roughness = (213.5*Re^-0.2841 + 1.082)/1000;
    Cf_Htail(i) = Htail_roughness;
end

M_0 = 0.5;
Z = (((2-M_0^2)*cosd(sweep_Htail))/sqrt(1-M_0^2*(cosd(sweep_Htail))^2))/10;
K1 = 1 + Z*tc_Htail + 100*tc_Htail;

k_Htail = K1;

Swet_Htail = 1.02*2*S_ref_Htail;

% Pylons
for i = 1:1:i_max
    v = v_min + (i-1)*v_step;
    Re = (rho*v*MAC_pylons)/u;
    pylons_roughness = (213.5*Re^-0.2841 + 1.082)/1000;
    Cf_pylons(i) = pylons_roughness;
end

M_0 = 0.5;
Z = ((2-M_0^2)*cosd(sweep_pylons))/sqrt(1-M_0^2*(cosd(sweep_pylons))^2);
K1 = 1 + Z*tc_pylons + 100*tc_pylons;

k_pylons = K1;

Swet_pylons = 1.02*2*S_ref_pylons;

% Nacelles
for i = 1:1:i_max
    v = v_min + (i-1)*v_step;
    Re = (rho*v*l_n)/u;
    nacelles_roughness = (213.5*Re^-0.2841 + 1.082)/1000;
    Cf_nacelles(i) = nacelles_roughness;
end

K2 = 5.3933*exp(-0.9111*ratio_n) + 1.3965*exp(-0.0279*ratio_n);
k_n = K2;

S_wet_n = S_wet_n;

%% Parasite Drag
for i = 1:1:i_max
    v1 = v_min + (i-1)*v_step;
    Cd_p1 = ((k_wing*Cf_wing(i)*Swet_wing) + (k_f*Cf_f(i)*Swet_f) + (k_Vtail*Cf_Vtail(i)*Swet_Vtail) + (k_Htail*Cf_Htail(i)*Swet_Htail) + (k_pylons*Cf_pylons(i)*Swet_pylons) + (k_n*Cf_nacelles(i)*S_wet_n))/(10*Sref_wing);
    Cd_p(i) = Cd_p1;
    D_p1 = 1/2*rho*v1^2*Cd_p1*Sref_wing;
    D_p(i) = D_p1;
    v(i) = v1;
end

figure(4)
plot(v,D_p)
title('Parasite Drag vs. Velocity')
xlabel('Velocity (ft/s)')
ylabel('Drag (lbs)')
grid on

%% Induced Drag

% Efficiency Digitization
AR = AR;
L = W;
load('efficiency_for_real.mat')



for i = 1:1:i_max
    v1 = v_min + (i-1)*v_step;
    C_l = (2*L)/(rho*v1^2*Sref_wing);
    e = efficiency_brand_new(Cd_p(i));
    D_i1 = ((C_l^2)/(pi*AR*e))*100000;
    D_i(i) = D_i1;
end

figure(5)
plot(v,D_i)
title('Induced Drag vs. Velocity')
ylabel('D_i (lbs)')
xlabel('Velocity (ft/s)')

%% Total Drag and L/D

D = D_i + D_p;

figure(6)
plot(v,D_p)
hold on
plot(v,D_i,'-r')
plot(v,D)
hold off
title('Total, Induced, and Parasitic Drag vs. Velocity')
xlabel('Velocity (ft/s)')
ylabel('Drag (lbs)')
legend('Parasitic Drag','Induced Drag','Total Drag')
grid on

L_D = L./D;

figure(7)
plot(v,L_D)
title('L/D ratio vs. Velocity')
xlabel('Velocity (ft/s)')
ylabel('L/D')












