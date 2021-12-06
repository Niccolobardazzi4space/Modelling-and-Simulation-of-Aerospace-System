% Modelling and Simulation of Aerospace Systems (2021/2022)
% Assignment # 2
% Author: Niccolò Bardazzi 10800456

%% EX 1
clearvars; close all; clc; 

data.J1 = 0.2;
data.J2 = 0.1;
data.T0 = 0.1;
data.tf = 10;

data.m = 20.e-3;
data.r = 5.e-5;
data.L = 28.e-3;
data.Ixx = 1/2*data.m*data.r^2;
data.G = 25.5e9;
data.k = data.G*data.Ixx/data.L;
data.b = 4;

data.theta0 = zeros(4,1);

data.h = 1/100;
N = (data.tf-0)/data.h;
data.t = linspace(0, data.tf, N+1);
data.options_ode = odeset('reltol', 1e-6, 'abstol', 1e-6);
[t, y] = ode113(@(t,theta) reaction_wheel(t, theta, data), data.t, data.theta0, data.options_ode);

figure()
semilogy(t,y(:,1),'LineStyle','--','DisplayName','\theta_1','LineWidth',1.2), hold on
semilogy(t,y(:,2),'LineStyle','-.','DisplayName','\theta_2','LineWidth',1.2)
legend show, legend('FontSize',13)
title('Angle response',Interpreter='latex')
xlabel('t [s]', 'FontSize', 13,Interpreter='latex')
ylabel('$\theta [rad]$', 'FontSize', 13,Interpreter='latex')
grid on
ax = gca;
ax.FontSize = 13;

table = readtable('samples.txt');
measures = table2array(table(:,2:3));
time = table2array(table(:,1));
% Display of data in sample.txt:
figure(3)
plot(time, measures(:,1), time, measures(:,2))
grid on
title('Acceleration sampled from gyros',Interpreter='latex',FontSize=14)
xlabel('t [s]', 'FontSize', 13,Interpreter='latex')
ylabel('$\ddot{\theta} [rad/s^2]$', 'FontSize', 13,Interpreter='latex')
grid on
ax = gca;
ax.FontSize = 13;

x0 = [data.b data.k];
options = optimset('Display', 'iter');
x = fminsearch(@(x) obj(x, measures, data), x0, options);

figure()
plot(t,y(:,3),'LineStyle','--','DisplayName','$\dot{\theta_1}$','LineWidth',1.2), hold on
plot(t,y(:,4),'LineStyle','-.','DisplayName','$\dot{\theta_2}$','LineWidth',1.2)
legend show, legend('FontSize',13,Interpreter='latex')
title('Velocity response',Interpreter='latex')
xlabel('t [s]', 'FontSize', 13,Interpreter='latex')
ylabel('$\dot{\theta}$ [rad/s]', 'FontSize', 13,Interpreter='latex')
grid on
ax = gca;
ax.FontSize = 13;


%% EX 2
clearvars; close all; clc; 

% Accumulator
sys.rho_skydrol = 890;
% Adiabatic transformation:
sys.gamma = 1.2;
sys.P0 = 21.e6;
sys.P_inf = 2.5e6;
sys.V_N_inf = 10.e-3;

% Section 23
sys.kA = 1.12;
sys.kcv = 2;
sys.D23 = 18.e-3; 
% assuming DA = D23
sys.DA = sys.D23;
sys.L23 = 2;
sys.f23 = 0.032;

% Actuator
sys.Dc = 50.e-3;
sys.Ac = pi*sys.Dc^2/4;
sys.Dr = 22.e-3;
sys.Ar = pi*sys.Dr^2/4;
sys.xmax = 200.e-3;
sys.m_piston = 2;
sys.k_piston = 120.e3;
sys.F0 = 1.e3;

% Tank
sys.P_T = 0.1e6;
sys.V_T0 = 1.e-3;

% Section 67
sys.D67 = sys.D23;
sys.L67 = 15;
sys.f67 = 0.035;
sys.k67 = 1.12;

% distributor
sys.dd = 5.e-3;
sys.kd = 12;

sys.t0 = 0;
sys.t_start = 1;
sys.t_close = 1.5;
sys.tf = 3;

sys.V_N_0 = sys.V_N_inf*(sys.P_inf/sys.P0); % isothermal transformation
sys.Vacc0 = sys.V_N_inf-sys.V_N_0;
x0 = [0; 0; sys.Vacc0; sys.V_T0];

% Variables of integration: - x(1): position of piston, x0(1) = 0;
%                           - x(2): velocity of piston, x0(2) = 0;
%                           - x(3): Vacc, Vacc0 = 10dm^3V_N_0;
%                           - x(4): VT, VT0;

options_ode = odeset('reltol', 1e-8, 'abstol', 1e-8, 'event', @max_stroke);
[t1, y1, te, ye, ie] = ode15s(@hydraulic_sys, [sys.t0 sys.tf], x0, options_ode, sys);
ye(2) = 0;
options_ode = odeset('reltol', 1e-8, 'abstol', 1e-8);
[t2, y2] = ode15s(@(t,y) hydraulic_sys(t, y, sys), [te sys.tf], ye, options_ode);
y = [y1; y2];
t = [t1; t2];
fprintf('Time to reach the maximum stroke: %.3f sec \n', te)

P = zeros(length(y),8);
dxdt = zeros(length(y),4);
for i = 1:length(y)
    [dxdt(i,:),P(i,1),P(i,2:end)] = hydraulic_sys(t(i), y(i,:), sys);
end

figure()
plot(t,y(:,1),'LineStyle','-','DisplayName','position','LineWidth',1.2), hold on
plot(t,y(:,2),'LineStyle','--','DisplayName','velocity','LineWidth',1.2), hold on
legend show, legend('FontSize',13)
title('State of the piston')
xlabel('t [s]', 'FontSize', 13)
yyaxis left
ylabel('x [m]', 'FontSize', 13,'Color', [0 0.4470 0.7410])
yyaxis right
ylabel('v [m/s]', 'FontSize', 13,'Color', [0.8500 0.3250 0.0980])
ax = gca;
ax.YAxis(1).Color = [0 0.4470 0.7410];
ax.YAxis(2).Color = [0.8500 0.3250 0.0980];
grid on
ax.FontSize = 13;

figure()
plot(t,y(:,3),'LineStyle','-','DisplayName','V_{acc}','LineWidth',1.2), hold on
plot(t,y(:,4),'LineStyle','--','DisplayName','V_{tank}','LineWidth',1.2), hold on
legend show, legend('FontSize',13)
title('Volumes')
xlabel('t [s]', 'FontSize', 13)
ylabel('V [m^3]', 'FontSize', 13)
grid on
ax = gca;
ax.FontSize = 13;

figure()
plot(t,P(:,1),'LineStyle','-','DisplayName','p_{acc}','LineWidth',1.2), hold on
plot(t,P(:,2),'LineStyle','--','DisplayName','p_1','LineWidth',1.2), hold on
plot(t,P(:,3),'LineStyle','--','DisplayName','p_2','LineWidth',1.2), hold on
plot(t,P(:,4),'LineStyle','--','DisplayName','p_3','LineWidth',1.2), hold on
plot(t,P(:,5),'LineStyle','--','DisplayName','p_4','LineWidth',1.2), hold on
plot(t,P(:,6),'LineStyle','--','DisplayName','p_5','LineWidth',1.2), hold on
plot(t,P(:,7),'LineStyle','--','DisplayName','p_6','LineWidth',1.2), hold on
plot(t,P(:,8),'LineStyle','--','DisplayName','p_7','LineWidth',1.2), hold on
legend show, legend('FontSize',13)
title('Pressures')
xlabel('t [s]', 'FontSize', 13)
ylabel('p [Pa]', 'FontSize', 13)
grid on
ax = gca;
ax.FontSize = 13;


%% EX 3
clearvars; close all; clc; 

data.R1 = 1000;
data.R2 = 100;
data.L = 1.e-3;
data.C = 1.e-3;
data.Vc0 = [1; 0];
data.f = 5;

R2 = data.R2;
R1 = data.R1;
L = data.L;
C = data.C;
f = data.f;

syms x  

vg(x) = 0*x;
tf = 5;

options_ode = odeset('reltol', 1e-6, 'abstol', 1e-6);
data.voltage_source = 0;
[t, y] = ode15s(@(t,Vc) circuit_integration(t, Vc, data), [0 tf], data.Vc0, options_ode);
figure()
plot(t,y(:,1),'LineStyle','-','DisplayName','V_{C}'), hold on
grid on
title('Capacitor discharge',Interpreter='latex')
xlabel('t [s]',FontSize=12)
ylabel('V_c [V]',FontSize=12)
ax = gca;
ax.FontSize = 12; 

vg(x) = sin(2*pi*f*x)*atan(x);
dvgdt = diff(vg,x);
data.voltage_source = 1;
[t, y] = ode15s(@(t,Vc) circuit_integration(t, Vc, data), [0 tf], data.Vc0, options_ode);
[dVcdt, A] = circuit_integration(t(1), y(1,:)', data);
eigenvalues = eig(A);
fprintf('The eigenvalues of the system are:\n [ %.2f \n %.2f ]\n', eigenvalues(1),eigenvalues(2))
figure()
plot(t,y(:,1),'LineStyle','-','DisplayName','V_{C}','Color',[0.8500, 0.3250, 0.0980]), hold on
grid on
title('Voltage source: V(t) = $\sin(2\pi ft)\arctan(t)$ and capacitor discharge',Interpreter='latex')
xlabel('t [s]',FontSize=12)
ylabel('V_c [V]',FontSize=12)
ax = gca;
ax.FontSize = 12; 


%% EX 4 
clearvars; close all; clc; 

data.k1 = 200;       
data.k2 = 16.3;      
data.k3 = 0.04;      
data.k4 = 0.4;    
data.k5 = 45;        

data.c2 = 680;       data.rho2 = 1600;
data.c4 = 1900;       data.rho4 = 200;

data.l(1) = 7.e-3; % mm
data.l(2) = 20.e-3;
data.l(3) = 5.e-4;
data.l(4) = 25.e-3;
data.l(5) = 10.e-3;

data.D = 525.e-3;
data.L = 500.e-3;

data.T_start = 20;
data.T_fire = 1000;
data.tstart = 0;
data.ti = 1;
data.tf = 60;
data.tspan = [data.ti, data.tf];  N = 500;

% FINITE DIFFERENCES, CYLINDRICAL RESISTANCES
% 1 point for insulator and conductor layers
data.points = 1;
[T_1,t_1] = temperature_integration(data, N);
figure()
plot(t_1,T_1(:,1),'DisplayName','Inner','LineStyle','-','LineWidth',1.2), hold on
plot(t_1,T_1(:,2),'DisplayName','Layer 1','LineStyle','--','LineWidth',1.2), hold on
plot(t_1,T_1(:,3),'DisplayName','Layer 2','LineStyle','-.','LineWidth',1.2), hold on
plot(t_1,T_1(:,4),'DisplayName','Layer 3','LineStyle',':','LineWidth',1.2), hold on
plot(t_1,T_1(:,5),'DisplayName','Layer 4','LineStyle','-','LineWidth',1.2), hold on
plot(t_1,T_1(:,6),'DisplayName','Layer 5','LineStyle','--','LineWidth',1.2), hold on
plot(t_1,T_1(:,7),'DisplayName','Outer','LineStyle','-.','LineWidth',1.2), hold on
title('FD, Cylidrical resistances, 1 layer')
legend show, legend(FontSize=13)
xlabel('t [s]', FontSize=13)
ylabel('T [°C]', FontSize=13)
xlim([min(t_1) max(t_1)])
grid on
ax = gca;
ax.FontSize = 12; 

% 2 points for insulator and conductor layers
data.points = 2;
[T_2,t_2] = temperature_integration(data, N);
figure()
plot(t_2,T_2(:,1),'DisplayName','Inner','LineStyle','-','LineWidth',1.2), hold on
plot(t_2,T_2(:,2),'DisplayName','Layer 1','LineStyle','--','LineWidth',1.2), hold on
plot(t_2,T_2(:,3),'DisplayName','Layer 2, 1^{st} point','LineStyle','-.','LineWidth',1.2,'Color',[0.9290, 0.6940, 0.1250]), hold on
plot(t_2,T_2(:,4),'DisplayName','Layer 2, 2^{nd} point','LineStyle',':','LineWidth',1.2,'Color',[0.9290, 0.6940, 0.1250]), hold on
plot(t_2,T_2(:,5),'DisplayName','Layer 3','LineStyle','-','LineWidth',1.2,'Color',[0.4940, 0.1840, 0.5560]), hold on
plot(t_2,T_2(:,6),'DisplayName','Layer 4, 1^{st} point','LineStyle','--','LineWidth',1.2,'Color',[0.4660, 0.6740, 0.1880]), hold on
plot(t_2,T_2(:,7),'DisplayName','Layer 4, 2^{nd} point','LineStyle','-.','LineWidth',1.2,'Color',[0.4660, 0.6740, 0.1880]), hold on
plot(t_2,T_2(:,8),'DisplayName','Layer 5','LineStyle',':','LineWidth',1.2,'Color',[0.3010, 0.7450, 0.9330]), hold on
plot(t_2,T_2(:,9),'DisplayName','Outer','LineStyle','-','LineWidth',1.2,'Color',[0.6350, 0.0780, 0.1840]), hold on
title('FD, Cylidrical resistances, 2 layers')
legend show, legend(FontSize=13)
xlabel('t [s]', FontSize=13)
ylabel('T [°C]', FontSize=13)
xlim([min(t_2) max(t_2)])
ylim([0 1200])
grid on
ax = gca;
ax.FontSize = 12;

% INTEGRATION WITH ODE, PLANAR RESISTANCES
% 1 point for insulator and conductor layers
T0 = ones(7,1)*20;
data.R1 = data.l(1)/(2*data.k1*data.L*pi*data.D);
data.R2 = data.l(2)/(2*data.k2*data.L*pi*data.D);
data.R3 = data.l(3)/(2*data.k3*data.L*pi*data.D);
data.R4 = data.l(4)/(2*data.k4*data.L*pi*data.D);
data.R5 = data.l(5)/(2*data.k5*data.L*pi*data.D);
m_A2 = data.rho2*pi*data.L*data.D;
m_A4 = data.rho4*pi*data.L*data.D;
C(1) = m_A2*data.l(2)*data.c2;
C(2) = m_A4*data.l(4)*data.c4;
A = [                 0                                                 0                                                         0                                                       0
         1/((2*data.R1+data.R2)*C(1))        -1/((2*data.R1+data.R2)*C(1))-1/((2*data.R3+data.R2+data.R4)*C(1))   -1/((2*data.R3+data.R2+data.R4)*C(1))                                   0
                      0                                 1/((2*data.R3+data.R2+data.R4)*C(2))    -1/((2*data.R3+data.R2+data.R4)*C(2))-1/((2*data.R5+data.R4)*C(2))        -1/((2*data.R5+data.R4)*C(2)) 
                      0                                                 0                                                         0                                                       0               ];
eigenvalues = eig(A);
fprintf('The eigenvalues of the system are:\n [ %.4f \n %.4f \n %.4f \n %.4f]\n', eigenvalues(1),eigenvalues(2),eigenvalues(3),eigenvalues(4))

on = 1; off = 0;
data.ramp = on;
options_ode = odeset('reltol', 1e-6, 'abstol', 1e-6);
[t_1_ramp, T_1_ramp] = ode113(@(t,T) temperature_1(t, T, data), [data.tstart data.ti], T0, options_ode);
data.ramp = off;
[t_1_const, T_1_const] = ode113(@(t,T) temperature_1(t, T, data), [data.ti data.tf], T_1_ramp(end,:)', options_ode);
t_1 = [t_1_ramp; t_1_const];
T_1 = [T_1_ramp; T_1_const];
T_1(:,2) = (T_1(:,1)/data.R1+T_1(:,3)/(data.R1+data.R2))*((data.R1+data.R2)*data.R1)/(2*data.R1+data.R2);
T_1(:,4) = (T_1(:,5)/(data.R4+data.R3)+T_1(:,3)/(data.R2+data.R3))*(data.R4+data.R3)*(data.R2+data.R3)/(data.R2+2*data.R3+data.R4);
T_1(:,6) = (T_1(:,7)/data.R5+T_1(:,5)/(data.R4+data.R5))*(data.R5*(data.R4+data.R5))/(2*data.R5+data.R4);
figure()
plot(t_1,T_1(:,1),'DisplayName','Inner','LineStyle','-','LineWidth',1.2), hold on
plot(t_1,T_1(:,2),'DisplayName','Layer 1','LineStyle','--','LineWidth',1.2), hold on
plot(t_1,T_1(:,3),'DisplayName','Layer 2','LineStyle','-.','LineWidth',1.2), hold on
plot(t_1,T_1(:,4),'DisplayName','Layer 3','LineStyle',':','LineWidth',1.2), hold on
plot(t_1,T_1(:,5),'DisplayName','Layer 4','LineStyle','-','LineWidth',1.2), hold on
plot(t_1,T_1(:,6),'DisplayName','Layer 5','LineStyle','--','LineWidth',1.2), hold on
plot(t_1,T_1(:,7),'DisplayName','Outer','LineStyle','-.','LineWidth',1.2), hold on
title('ODE, planar resistances, 1 layer')
legend show, legend(FontSize=13)
xlabel('t [s]', FontSize=13)
ylabel('T [°C]', FontSize=13)
xlim([min(t_1) max(t_1)])
grid on
ax = gca;
ax.FontSize = 12; 

% 2 points for insulator and conductor layers
T0 = ones(9,1)*20;
data.ramp = on;
options_ode = odeset('reltol', 1e-6, 'abstol', 1e-6);
[t_2_ramp, T_2_ramp] = ode113(@(t,T) temperature_2(t, T, data), [data.tstart data.ti], T0, options_ode);
data.ramp = off;
[t_2_const, T_2_const] = ode113(@(t,T) temperature_2(t, T, data), [data.ti data.tf], T_2_ramp(end,:)', options_ode);
t_2 = [t_2_ramp; t_2_const];
T_2 = [T_2_ramp; T_2_const];
data.R2 = data.l(2)/(3*data.k2*data.L*pi*data.D);
data.R4 = data.l(4)/(3*data.k4*data.L*pi*data.D);
m_A2 = data.rho2*pi*data.L*data.D;
m_A4 = data.rho4*pi*data.L*data.D;
C(1) = m_A2*data.l(2)*data.c2/2;
C(2) = m_A4*data.l(4)*data.c4/2;
A = [                    0                                       0                                                    0                                                 0                                                   0                                                       0
     -1/((2*data.R1+data.R2)*C(1))              -1/((2*data.R1+data.R2)*C(1))-1/(data.R2*C(1))                 -1/(data.R2*C(1))                                        0                                                   0                                                       0
                         0                               -1/(data.R2*C(1))                          1/(data.R2*C(1))-1/((2*data.R3+data.R2+data.R4)*C(1))      -1/((2*data.R3+data.R2+data.R4)*C(1))                        0                                                       0
                         0                                       0                                       -1/((2*data.R3+data.R2+data.R4)*C(2))    -1/((2*data.R3+data.R2+data.R4)*C(2))-1/(data.R4*C(2))           -1/(data.R4*C(2))                                                0    
                         0                                       0                                                    0                                           -1/(data.R4*C(2))                            -1/(data.R4*C(2))-1/((2*data.R5+data.R4)*C(2))         -1/((2*data.R5+data.R4)*C(2))  
                         0                                       0                                                    0                                                 0                                                   0                                                       0         ];
eigenvalues = eig(A);
fprintf('The eigenvalues of the system are:\n [ %.4f \n %.4f \n %.4f \n %.4f \n %.4f \n %.4f]\n', eigenvalues(1),eigenvalues(2),eigenvalues(3),eigenvalues(4),eigenvalues(5),eigenvalues(6))

T_2(:,2) = (T_2(:,1)/data.R1+T_2(:,3)/(data.R1+data.R2))*((data.R1+data.R2)*data.R1)/(2*data.R1+data.R2);
T_2(:,5) = (T_2(:,6)/(data.R4+data.R3)+T_2(:,4)/(data.R2+data.R3))*(data.R4+data.R3)*(data.R2+data.R3)/(data.R2+2*data.R3+data.R4);
T_2(:,8) = (T_2(:,9)/data.R5+T_2(:,7)/(data.R4+data.R5))*(data.R5*(data.R4+data.R5))/(2*data.R5+data.R4);
figure()
plot(t_2,T_2(:,1),'DisplayName','Inner','LineStyle','-','LineWidth',1.2), hold on
plot(t_2,T_2(:,2),'DisplayName','Layer 1','LineStyle','--','LineWidth',1.2), hold on
plot(t_2,T_2(:,3),'DisplayName','Layer 2, 1^{st} point','LineStyle','-.','LineWidth',1.2,'Color',[0.9290, 0.6940, 0.1250]), hold on
plot(t_2,T_2(:,4),'DisplayName','Layer 2, 2^{nd} point','LineStyle',':','LineWidth',1.2,'Color',[0.9290, 0.6940, 0.1250]), hold on
plot(t_2,T_2(:,5),'DisplayName','Layer 3','LineStyle','-','LineWidth',1.2,'Color',[0.4940, 0.1840, 0.5560]), hold on
plot(t_2,T_2(:,6),'DisplayName','Layer 4, 1^{st} point','LineStyle','--','LineWidth',1.2,'Color',[0.4660, 0.6740, 0.1880]), hold on
plot(t_2,T_2(:,7),'DisplayName','Layer 4, 2^{nd} point','LineStyle','-.','LineWidth',1.2,'Color',[0.4660, 0.6740, 0.1880]), hold on
plot(t_2,T_2(:,8),'DisplayName','Layer 5','LineStyle',':','LineWidth',1.2,'Color',[0.3010, 0.7450, 0.9330]), hold on
plot(t_2,T_2(:,9),'DisplayName','Outer','LineStyle','-','LineWidth',1.2,'Color',[0.6350, 0.0780, 0.1840]), hold on
title('ODE, planar resistances, 2 layers')
legend show, legend(FontSize=13)
xlabel('t [s]', FontSize=13)
ylabel('T [°C]', FontSize=13)
grid on
ylim([0 1200])
ax = gca;
ax.FontSize = 12;


%% functions

% EX 1
function dthetadt = reaction_wheel(~, theta, data)
%
%     reaction_wheel - Integration of the mechanical system (reaction 
%                      wheel) in Ex.1
%
%     DESCRIPTION:
%       This function gives as output the derivative vector to be
%       integrated of the model of Ex.1
%
%     PROTOTYPE:
%        dthetadt = reaction_wheel(~, theta, data)
%          
%     INPUT:
%       t                  Integrator's time
%       theta [4,1]        Integrated vector at each time instant
%       data               Struct with fields: 
%                          - J1, Inertia of the first disk
%                          - J2, Inertia of the second disk
%                          - k, elastic constant of the spring
%                          - b, quadratic friction coefficient
%                          - T0, torque provided as input to the system
%     
%     OUTPUT:
%       dthetadt [4,1]     Vector to be integrated by ODE schemes
%     
%     CALLED FUNCTIONS:
%      -
% 
%     LAST UPDATED:
%      18/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

dthetadt = zeros(4,1);
dthetadt(1:2) = theta(3:4);
dthetadt(3) = -data.k*(theta(1)-theta(2))/data.J1;
dthetadt(4) = (-data.k*(theta(2)-theta(1))-sign(theta(4))*data.b*theta(4)^2+data.T0)/data.J2;
end

function fx = obj(x, measures, data)
%
%     obj - Objective function to be minimized
%
%     DESCRIPTION:
%       This function gives as output the objective function to be 
%       minimized in Ex.1 to get the correct values of k, b corresponding
%       to the accelerations sampled in sample.txt file
%
%     PROTOTYPE:
%        fx = obj(x, measures, data)
%          
%     INPUT:
%       x [2,1]            Minimizing variable
%       measures [n,2]     Measurements from sample.txt file
%       data               Struct with fields: 
%                          - J1, Inertia of the first disk
%                          - J2, Inertia of the second disk
%                          - k, elastic constant of the spring
%                          - b, quadratic friction coefficient
%                          - T0, torque provided as input to the system
%     
%     OUTPUT:
%       fx                 Scalar objective of the minimization
%     
%     CALLED FUNCTIONS:
%      -
% 
%     LAST UPDATED:
%      18/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

data.b = x(1);
data.k = x(2);
[t, theta] = ode15s(@(t,theta) reaction_wheel(t, theta, data), data.t, data.theta0, data.options_ode);
y(:,1) = -data.k*(theta(:,1)-theta(:,2))/data.J1;
y(:,2) = (-data.k*(theta(:,2)-theta(:,1))-sign(theta(:,4)).*data.b.*theta(:,4).^2+data.T0)/data.J2;
figure(2)
plot(t, y(:,1), t, y(:,2), t, measures(:,1), t, measures(:,2))
grid on
title('Acceleration with $k$ and $b$ found','Interpreter','latex',FontSize=14)
xlabel('t [s]', 'FontSize', 13,Interpreter='latex')
ylabel('$\ddot{\theta} [rad/s^2]$', 'FontSize', 13,Interpreter='latex')
grid on
ax = gca;
ax.FontSize = 13;
fx = norm([norm((y(:,1)-measures(:,1))); norm((y(:,2)-measures(:,2)))]);
end

% EX 2
function [value, isterminal, direction] = max_stroke(~, y, sys)
%
%     max_stroke - Event of Ex.2 when the piston reach the maximum stroke
%
%     DESCRIPTION:
%       This function gives as output variables needed from the ODE option
%       event to stop the integration whenever the piston reaches the
%       maximum stroke
%
%     PROTOTYPE:
%        [value, isterminal, direction] = max_stroke(~, y, sys)
%          
%     INPUT:
%       t                  Integrator's time
%       y [4,1]            Integrated vector at each time instant
%       sys                Struct with fields: 
%                          - xmax, maximum stroke 
%     
%     OUTPUT:
%       value              Stopping criteria value, x == x_max
%       isterminal         Stop integration when event occurs, 1
%       direction          Generic direction for stopping criteria, 0
%     
%     CALLED FUNCTIONS:
%      -
% 
%     LAST UPDATED:
%      18/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

value = sys.xmax-y(1);
isterminal = 1;
direction = 0;
end

function [dxdt,PA,P] = hydraulic_sys(t, x, sys)
%
%     hydraulic_sys - Mathematical model of the hydraulic system of Ex.2
%
%     DESCRIPTION:
%       This function gives as output the derivative vector to be
%       integrated, the pressure of the accumulator and the pressures in
%       the whole system
%
%     PROTOTYPE:
%        [dxdt,PA,P] = hydraulic_sys(t, x, sys)
%          
%     INPUT:
%       t                   Integrator's time
%       x [4,1]             Integrated vector at each time instant
%       sys                 Struct with fields: 
%                           - sys.rho_skydrol, fluid density
%                           - sys.gamma, adiabatic transformation coeff.
%                           - sys.P0, final isothermal pressure
%                           - sys.V_N_0, final isothermal volume
%                           - sys.P_inf, initial isothermal pressure
%                           - sys.V_N_inf, initial isothermal volume
%                           - sys.Vacc0, initial accumulation volume
%                           - sys.kA, coeff. of conc. pressure loss in A
%                           - sys.kcv, coeff. of conc. pressure loss in 72
%                           - sys.D23, diameter of pipe in 23 section
%                           - sys.DA, diameter of pipe in A section
%                           - sys.L23, length of pipe in 23 section
%                           - sys.f23, coeff. of distr. pressure loss in 23
%                           - sys.Dc, diameter of the cylinder in actuator
%                           - sys.Ac, area of the cylinder in actuator
%                           - sys.Dr, diameter of the rod in actuator
%                           - sys.Ar, area of the rod in actuator
%                           - sys.xmax, maximum stroke
%                           - sys.m_piston, mass of the piston
%                           - sys.k_piston, spring of the actuator
%                           - sys.F0, spring pre-load
%                           - sys.P_T, tank pressure
%                           - sys.V_T0, tank volume at initial time
%                           - sys.D67, diameter of pipe in 67 section
%                           - sys.L67, length of pipe in 67 section
%                           - sys.f67, coeff. of distr. pressure loss in 67
%                           - sys.k67, coeff. of conc. pressure loss in 67
%                           - sys.dd, diameter of the distributor 
%                           - sys.kd, coeff. of conc. pressure loss in
%                             ditributor
%                           - sys.t0, integration time starting
%                           - sys.t_start, starting time for valve opening
%                           - sys.t_close, ending time for valve closure
%                           - sys.tf, integration time ending
%     
%     OUTPUT:
%       dxdt [4,1]         Vector to be integrated by ODE schemes
%       PA                 Accumulator pressure
%       P [7,1]            Vector of pressures along circuit of Ex.2
%     
%     CALLED FUNCTIONS:
%      -
% 
%     LAST UPDATED:
%      18/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

    z = ramp(t, sys);
    % Different type of valve:
%     alpha = 2*acos(1-abs(z));
%     A_valve = sys.dd^2/4*(alpha-sin(alpha));
    z = z*pi/2;
    alpha = 2*acos(1-cos(z));
    A_valve = sys.dd^2/4*pi-(sys.dd^2*alpha/8-sys.dd^2*sin(alpha/2)*cos(alpha/2)/4)*2;
    dxdt = zeros(4,1);
    dxdt(1) = x(2);
    V_N_0 = sys.V_N_inf-x(3);
    % Accumulator
    PA = sys.P0*(sys.V_N_0/V_N_0)^sys.gamma; % adiabatic transformation
    % losses
    dxdt(3) = -dxdt(1)*(sys.Ac-sys.Ar);
    Qacc = dxdt(3);
    vA = Qacc/(pi*sys.DA^2/4);
    dp1 = 1/2*sys.rho_skydrol*sys.kA*vA^2;
    dp2 = 1/2*sys.rho_skydrol*sys.kcv*vA^2;
    P(1) = PA-dp1;
    P(2) = P(1)-dp2;
    v23 = Qacc/(pi*sys.D23^2/4);
    dp23 = sys.f23*sys.L23*sys.rho_skydrol*v23^2/(2*sys.D23);
    P(3) = P(2)-dp23;

    if z == 0 || A_valve<1.e-9
        P(1:3) = sys.P0; 
        P(4) = (sys.P_T*(sys.Ac-sys.Ar)+sys.F0)/sys.Ac;
        PA = sys.P0;
        P(6:7) = sys.P_T;
        P(5) = sys.P_T;
    end
    
    vd1 = Qacc/A_valve;

    if A_valve <1.e-9
        vd1 = 0;
    end
    dp_distributor1 = 1/2*sys.rho_skydrol*sys.kd*vd1^2;
    P(4) = P(3)-dp_distributor1;
    
    % Lower part
    dxdt(4) = dxdt(1)*(sys.Ac-sys.Ar);
    QT = dxdt(4);
    v67 = QT/(pi*sys.D67^2/4);
    P(7) = sys.P_T+1/2*sys.k67*v67^2;
    dp67 = sys.f67*sys.L67*sys.rho_skydrol*v67^2/(2*sys.D67);
    P(6) = P(7)+dp67;
    vd2 = QT/A_valve;

    if A_valve <1.e-9
        vd2 = 0;
    end

    dp_distributor2 = 1/2*sys.rho_skydrol*sys.kd*vd2^2;
    P(5) = P(6)+dp_distributor2;
    
    % solve piston to go back
    F = -sys.F0-sys.k_piston*x(1)+P(4)*sys.Ac-P(5)*(sys.Ac-sys.Ar);
    dxdt(2) = F/sys.m_piston;
    if x(1)>sys.xmax
        x(1) = sys.xmax;
        dxdt = zeros(4,1);
    end

    if x(1)<0
        dxdt = zeros(4,1);
    end

    if z == 0 || A_valve<1.e-9
        dxdt = zeros(4,1);
    end

end

function z = ramp(t, sys)
%
%     ramp - Linear valve opening
%
%     DESCRIPTION:
%       This function control the opening area of the valve in Ex.2
%
%     PROTOTYPE:
%        z = ramp(t, sys)
%          
%     INPUT:
%       t                  Integrator's time
%       sys                Struct with fields: 
%                          - sys.t0, integration time starting
%                          - sys.t_start, starting time for valve opening
%                          - sys.t_close, ending time for valve closure
%                          - sys.tf, integration time ending
%     
%     OUTPUT:
%       z                  Scalar between 0 and 1
%     
%     CALLED FUNCTIONS:
%      -
% 
%     LAST UPDATED:
%      18/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

if t<sys.t_start
    z = 0;
elseif t>sys.t_close
    z = 1;
else
    z = (t-sys.t_start)/(sys.t_close-sys.t_start);
end

end

% EX 3
function [dVcdt, A] = circuit_integration(t, Vc, data)
%
%     circuit_integration - Integration of the circuit in Ex.3
%
%     DESCRIPTION:
%       This function gives as output the derivative vector to be
%       integrated and the matrix whose eigenvalues determines the dynamic
%       behavior of the system of the model of Ex.3
%
%     PROTOTYPE:
%        [dVcdt, A] = circuit_integration(t, Vc, data)
%          
%     INPUT:
%       t                  Integrator's time
%       Vc [2,1]           Integrated vector at each time instant
%       data               Struct with fields: 
%                          - R1, value of R1 resistance
%                          - R2, value of R2 resistance
%                          - L, value of the inductance
%                          - C, value of the capacitance
%                          - f, freuency of the generator
%     
%     OUTPUT:
%       dVcdt [2,1]        Vector to be integrated by ODE schemes
%       A [2,2]            Matrix of dynamics
%     
%     CALLED FUNCTIONS:
%      -
% 
%     LAST UPDATED:
%      18/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

R2 = data.R2;
R1 = data.R1;
L = data.L;
C = data.C;
f = data.f;

dVcdt = zeros(2,1);
A = [           0                         1     
     -1/(L*C+R2*C*L/R1)  -(R2*C+L/R1)/(L*C+R2*C*L/R1)];

switch data.voltage_source
    case 1
        b = [                                                 0
            (L/R1*(sin(10*pi*t)/(t^2 + 1) + 10*pi*cos(10*pi*t)*atan(t))+sin(2*pi*f*t)*atan(t))/(L*C+R2*C*L/R1)];
        dVcdt = A*Vc+b;
    case 0
        dVcdt = A*Vc;
end
end

% EX 4
function [T,t] = temperature_integration(data, N)
%
%     temperature_integration - Integration of the system with Finite
%                               Differences
%
%     DESCRIPTION:
%       This function integrate the thermal system of Ex.4 performing a
%       with forward finite differences
%
%     PROTOTYPE:
%        [T,t] = temperature_integration(T_start, T_fire, data, N)
%          
%     INPUT:
%       data                Struct with fields: 
%                           - data.k1, conductivity of layer 1       
%                           - data.k2, conductivity of layer 2       
%                           - data.k3, conductivity of layer 3        
%                           - data.k4, conductivity of layer 4      
%                           - data.k5, conductivity of layer 5        
%                           - data.c2, capacity of layer 2       
%                           - data.rho2, density of layer 2
%                           - data.c4, capacity of layer 4     
%                           - data.rho4, density of layer 4
%                           - data.l, vector with the lenghts of all the
%                                     layers 
%                           - data.D, diameter of the rocket engine
%                           - data.L, length of the rocket engine
%                           - data.T_start, initial temperature of all
%                                           layers
%                           - data.T_fire, flame temperature inside rocket
%                           - data.tstart, firing starting time
%                           - data.ti, time when temperature inside rocket 
%                                      reaches T_fire                                      
%                           - data.tf, final integration time
%                           - data.points, # of point in capacitance layers
%       N                   # of time span intervals
%     
%     OUTPUT:
%       T                   Temperature vector whose length is
%                           (5+2*data.points,1)
%       t                   Integration time
%     
%     CALLED FUNCTIONS:
%       cond_row, cap_row
% 
%     LAST UPDATED:
%      18/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

k1 = data.k1;
k2 = data.k2;
k3 = data.k3;
k4 = data.k4;
k5 = data.k5;
tspan = data.tspan;
D = data.D;
L = data.L;
rho2 = data.rho2;
rho4 = data.rho4;
c2 = data.c2;
c4 = data.c4;
points = data.points;
T_start = data.T_start;
T_fire = data.T_fire;
t = linspace(0, data.tspan(2), N+2);

l = zeros(3+data.points*2,1);
for i = 1:length(data.l)
    if i == 2 || i == 4
        l(i) = data.l(i)/(points+1);
    else
        l(i) = data.l(i)/2;
    end
end

dt = (tspan(2)-0)/N;

A = zeros(length(l)+2,length(l)+2);
b = zeros(length(l)+2,1);
T = ones(N+1,length(l)+2)*T_start;
A(1,1) = 1;
b(1) = 0;
A(length(l)+2, length(l)+2) = 1;
b(end) = T_start;

if data.points == 1
    for i = 1:N+1
        % A matrix
        A(2,1:3) = cond_row([k1, k1, k2], L, D, [0, l(1), l(2)]);
        [A(3,2:4), b(3)] = cap_row([k1, k2, k3], dt, rho2, c2, [l(1) l(2) l(3)], D+l(1), L, T(i,3));
        A(4,3:5) = cond_row([k2, k3, k4], L, D+sum(l(1:2)), [l(2), l(3), l(4)]);
        [A(5,4:6), b(5)] = cap_row([k3, k4, k5], dt, rho4, c4, [l(3) l(4) l(5)], D+sum(l(1:3)), L, T(i,5));
        A(6,5:7) = cond_row([k4, k5, k5], L, D+sum(l(1:4)), [l(4) l(5) l(5)]);
        
        if t(i) >= data.tspan(1)
            b(1) = T_fire;
        else
            b(1) = 20+T_fire*t(i)/(data.tspan(1));
        end
        T(i+1,:) = A\b;
    end
else 
    for i = 1:N+1
        % A matrix
        A(2,1:3) = cond_row([k1, k1, k2], L, D+l(1), [0, l(1), l(2)]);
        [A(3,2:4), b(3)] = cap_row([k1, k2, k2], dt, rho2, c2, [l(1) l(2) l(2)], D+(l(1)*2+l(2)), L, T(i,3));
        [A(4,3:5), b(4)] = cap_row([k2, k2, k3], dt, rho2, c2, [l(2) l(2) l(3)], D+(l(1)*2+l(2)*2), L, T(i,4));
        A(5,4:6) = cond_row([k2, k3, k4], L, D+(l(1)*2+l(2)*3+l(3)), [l(2), l(3), l(4)]);
        [A(6,5:7), b(6)] = cap_row([k3, k4, k4], dt, rho4, c4, [l(3) l(4) l(4)], D+(l(1)*2+l(2)*3+l(3)*2+l(4)), L, T(i,6));
        [A(7,6:8), b(7)] = cap_row([k4, k4, k5], dt, rho4, c4, [l(4) l(4) l(5)], D+(l(1)*2+l(2)*3+l(3)*2+l(4)*2), L, T(i,7));
        A(8,7:9) = cond_row([k4, k5, k5], L, D+(l(1)*2+l(2)*3+l(3)*2+l(4)*3+l(5)), [l(4) l(5) l(5)]);          

        if t(i) >= data.tspan(1)
            b(1) = T_fire;
        else
            b(1) = 20+T_fire*t(i)/(data.tspan(1));
        end
        T(i+1,:) = A\b;
    end
end
end

function A = cond_row(k, L, D, thickness)
%
%     cond_row - conductive row for FD integration
%
%     DESCRIPTION:
%       This function gives as output a reistance vector to be placed in 
%       the row corresponding to a conductive layer inside the tridiagonal 
%       matrix of forward finite difference integration for a cylindrical
%       body thermal system
%
%     PROTOTYPE:
%        A = cond_row(k, L, D, thickness)
%          
%     INPUT:
%       k [3,1]             Conductivity vector of the previous, present
%                           and successive layer
%       L                   Length of the cylinder
%       D                   Diameter of the first layer
%       thickness [3,1]     Thickness vector of the previous, present and
%                           successive layer
%
%     OUTPUT:
%       A                   Reistance vector of the conductive layer
%     
%     CALLED FUNCTIONS:
%      -
% 
%     LAST UPDATED:
%      18/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

A = zeros(1,3);
R1 = log((D/2+thickness(1))./(D/2))./(2*pi*k(1)*L)+log((D/2+thickness(1)+thickness(2))./(D/2+thickness(1)))./(2*pi*k(2)*L);
R2 = log((D/2+thickness(1)+thickness(2)*2)./(D/2+thickness(1)+thickness(2)))./(2*pi*k(2)*L)+log((D/2+thickness(1)+thickness(2)*2+thickness(3))./(D/2+thickness(1)+thickness(2)*2))./(2*pi*k(3)*L);
A(2) = 1/R1+1/R2;
A(1) = -1/R1;
A(3) = -1/R2;
end

function [B, b] = cap_row(k, dt, rho, c, thickness, D, L, T) 
%
%     cap_row - conductive row for FD integration
%
%     DESCRIPTION:
%       This function gives as output a reistance vector to be placed in 
%       the row corresponding to a conductive layer inside the tridiagonal 
%       matrix of forward finite difference integration for a cylindrical
%       body thermal system
%
%     PROTOTYPE:
%        [B, b] = cap_row(k, dt, rho, c, thickness, D, L, T) 
%          
%     INPUT:
%       k [3,1]             Conductivity vector of the previous, present
%                           and successive layer
%       dt                  Time difference between previous and present
%                           time instant
%       rho                 Density of the capacitance layer
%       c                   Capacity of the layer
%       thickness [3,1]     Thickness vector of the previous, present and
%                           successive layer
%       D                   Diameter of the first layerù
%       L                   Length of the cylinder
%       T                   Temperature of the point at the previous time
%                           instant
%
%     OUTPUT:
%       B                   Reistance vector of the capacitive layer
%       b                   Temperature vector
%     
%     CALLED FUNCTIONS:
%      -
% 
%     LAST UPDATED:
%      18/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

B = zeros(1,3);
M = rho*pi*(((D+thickness(1)*2+thickness(2)*2)/2)^2-((D+thickness(1)*2)/2)^2);
R1 = log((D/2+thickness(1))./(D/2))./(2*pi*k(1)*L)+log((D/2+thickness(1)+thickness(2))./(D/2+thickness(1)))./(2*pi*k(2)*L);
R2 = log((D/2+thickness(1)+thickness(2)*2)./(D/2+thickness(1)+thickness(2)))./(2*pi*k(2)*L)+log((D/2+thickness(1)+thickness(2)*2+thickness(3))./(D/2+thickness(1)+thickness(2)*2))./(2*pi*k(3)*L);
B(2) = c*M*R1*R2/dt+R1+R2;                                                      
B(3) = -R1;
B(1) = -R2;
b = T*R1*R2*M*c/dt;
end

function dxdt = temperature_1(~, T, data)
%
%     temperature_1 - Integration of the circuit in Ex.4
%
%     DESCRIPTION:
%       This function gives as output the derivative vector to be
%       integrated by an ODE scheme representing the model of Ex.4 with 1
%       node per capacitance layers.
%
%     PROTOTYPE:
%        dxdt = temperature_1(~, T, data)
%          
%     INPUT:
%       t                  Integrator's time
%       T [7,1]            Integrated vector at each time instant
%       data                Struct with fields: 
%                           - data.R1, resistance of layer 1       
%                           - data.R2, resistance of layer 2       
%                           - data.R3, resistance of layer 3       
%                           - data.R4, resistance of layer 4    
%                           - data.R5, resistance of layer 5       
%                           - data.c2, capacity of layer 2       
%                           - data.rho2, density of layer 2
%                           - data.c4, capacity of layer 4     
%                           - data.rho4, density of layer 4
%                           - data.l, vector with the lenghts of all the
%                                     layers 
%                           - data.D, diameter of the rocket engine
%                           - data.L, length of the rocket engine
%                           - data.T_start, initial temperature of all
%                                           layers
%                           - data.T_fire, flame temperature inside rocket
%                           - data.tstart, firing starting time
%                           - data.ti, time when temperature inside rocket 
%                                      reaches T_fire                                      
%     
%     OUTPUT:
%       T                   Temperature vector whose length is
%                           (5+2*data.points,1)
%       t                   Integration time
%     
%     CALLED FUNCTIONS:
%       cond_row, cap_row
% 
%     LAST UPDATED:
%      18/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

m_A2 = data.rho2*pi*data.L*data.D;
m_A4 = data.rho4*pi*data.L*data.D;
C(1) = m_A2*data.l(2)*data.c2;
C(2) = m_A4*data.l(4)*data.c4;
dxdt = zeros(7,1);
if data.ramp == 1
    dxdt(1) = (data.T_fire-data.T_start)/(data.ti-data.tstart);
end
dxdt(3) = (T(1)-T(3))/((2*data.R1+data.R2)*C(1))-(T(3)-T(5))/((2*data.R3+data.R2+data.R4)*C(1));
dxdt(5) = (T(3)-T(5))/((2*data.R3+data.R2+data.R4)*C(2))-(T(5)-T(7))/((2*data.R5+data.R4)*C(2));
end

function dxdt = temperature_2(~, T, data)
%
%     temperature_2 - Integration of the circuit in Ex.4
%
%     DESCRIPTION:
%       This function gives as output the derivative vector to be
%       integrated by an ODE scheme representing the model of Ex.4 with 2
%       nodes per capacitance layers.
%
%     PROTOTYPE:
%        dxdt = temperature_2(~, T, data)
%          
%     INPUT:
%       t                  Integrator's time
%       T [9,1]            Integrated vector at each time instant
%       data                Struct with fields: 
%                           - data.R1, resistance of layer 1       
%                           - data.R2, resistance of layer 2       
%                           - data.R3, resistance of layer 3       
%                           - data.R4, resistance of layer 4       
%                           - data.R5, resistance of layer 5       
%                           - data.c2, capacity of layer 2       
%                           - data.rho2, density of layer 2
%                           - data.c4, capacity of layer 4     
%                           - data.rho4, density of layer 4
%                           - data.l, vector with the lenghts of all the
%                                     layers 
%                           - data.D, diameter of the rocket engine
%                           - data.L, length of the rocket engine
%                           - data.T_start, initial temperature of all
%                                           layers
%                           - data.T_fire, flame temperature inside rocket
%                           - data.tstart, firing starting time
%                           - data.ti, time when temperature inside rocket 
%                                      reaches T_fire                                      
%     
%     OUTPUT:
%       T                   Temperature vector whose length is
%                           (5+2*data.points,1)
%       t                   Integration time
%     
%     CALLED FUNCTIONS:
%       cond_row, cap_row
% 
%     LAST UPDATED:
%      18/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

m_A2 = data.rho2*pi*data.L*data.D;
m_A4 = data.rho4*pi*data.L*data.D;
C(1) = m_A2*data.l(2)*data.c2/2;
C(2) = m_A4*data.l(4)*data.c4/2;
dxdt = zeros(9,1);
if data.ramp == 1
    dxdt(1) = (data.T_fire-data.T_start)/(data.ti-data.tstart);
end
dxdt(3) = (T(1)-T(3))/((2*data.R1+data.R2)*C(1))-(T(3)-T(4))/(data.R2*C(1));
dxdt(4) = (T(3)-T(4))/(data.R2*C(1))-(T(4)-T(6))/((2*data.R3+data.R2+data.R4)*C(1));
dxdt(6) = (T(4)-T(6))/((2*data.R3+data.R2+data.R4)*C(2))-(T(6)-T(7))/(data.R4*C(2));
dxdt(7) = (T(6)-T(7))/(data.R4*C(2))-(T(7)-T(9))/((2*data.R5+data.R4)*C(2));
end


