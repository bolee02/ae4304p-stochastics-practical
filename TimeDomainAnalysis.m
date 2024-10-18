% Filename : TimeDomainAnalysis (originally examp81c.m)
% From lecture notes ae4-304
%
% Bo Lee (5225604)
%%

% RUN PREVIOUS SCRIPTS
StabilityAnalysis;
close all;

dt = 0.01; T  = 60; t = 0:dt:T; N = length(t);
nn = zeros(1,N);

% TURBULENCE INPUTS
u_g = randn(1,N)/sqrt(dt);    % sqrt(dt) because of lsim characteristics
w_g = randn(1,N)/sqrt(dt);

% INPUT VECTORS
u1 = [nn' nn' u_g' nn'  nn'];     
u2 = [nn' nn' nn'  w_g' nn'];

% DEFINE OUTPUT MATRICES
C = [1 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0 0;
     V*yb V*yphi V*yp V*(yr+2*V/b) 0 0 0 0 V*ybg 0];
   
D = [0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0];

% RESPONSE to u_g
y_u = lsim(A_cl,B,C,D,u1,t);

% RESPONSE to w_g
y_w = lsim(A_cl,B,C,D,u2,t);

% PLOT RESULTS for y_u
figure(4)
grid on
plot(t, y_u)
xlabel('Time [s]')
ylabel('Response - y_u')
set(gcf,'color','white')
legend({'\beta_u', '\phi_u', 'pb/2V_u', 'rb/2V_u', 'a_{y_u}'})
export_fig('timedomain_y_u', '-png', '-r300', '-nocrop')

figure(5)
grid on
plot(t, y_u(:,5))
xlabel('Time [s]')
ylabel('Lateral acceleration - a_{y_u} [m/s^2]')
set(gcf,'color','white')
legend('a_{y_u}')
export_fig('timedomain_ay_u', '-png', '-r300', '-nocrop')

figure(6)
grid on
plot(t, y_u(:,1))
xlabel('Time [s]')
ylabel('Sidelip - \beta_u [rad]')
set(gcf,'color','white')
legend('\beta_u')
export_fig('timedomain_beta_u', '-png', '-r300', '-nocrop')

figure(7)
grid on
plot(t, y_u(:,4))
xlabel('Time [s]')
ylabel('rb/2V_u [rad]')
set(gcf,'color','white')
legend('rb/2V_u')
export_fig('timedomain_rb_u', '-png', '-r300', '-nocrop')

figure(8)
grid on
plot(t, y_u(:,2))
xlabel('Time [s]')
ylabel('Roll - \phi_u [rad]')
set(gcf,'color','white')
legend('\phi_u')
export_fig('timedomain_phi_u', '-png', '-r300', '-nocrop')

figure(9)
grid on
plot(t, y_u(:,3))
xlabel('Time [s]')
ylabel('pb/2V_u [rad]')
set(gcf,'color','white')
legend('pb/2V_u')
export_fig('timedomain_pb_u', '-png', '-r300', '-nocrop')

% PLOT RESULTS for y_w
figure(10)
grid on
plot(t, y_w)
xlabel('Time [s]')
ylabel('Response - y_w')
set(gcf,'color','white')
legend({'\beta_w', '\phi_w', 'pb/2V_w', 'rb/2V_w', 'a_{y_w}'})
export_fig('timedomain_y_w', '-png', '-r300', '-nocrop')

figure(11)
grid on
plot(t, y_w(:,1))
xlabel('Time [s]')
ylabel('Sidelip - \beta_w [rad]')
set(gcf,'color','white')
legend('\beta_w')
export_fig('timedomain_beta_w', '-png', '-r300', '-nocrop')

figure(12)
grid on
plot(t, y_w(:,2))
xlabel('Time [s]')
ylabel('Roll - \phi_w [rad]')
set(gcf,'color','white')
legend('\phi_w')
export_fig('timedomain_phi_w', '-png', '-r300', '-nocrop')

figure(13)
grid on
plot(t, y_w(:,3))
xlabel('Time [s]')
ylabel('pb/2V_w [rad]')
set(gcf,'color','white')
legend('pb/2V_w')
export_fig('timedomain_pb_w', '-png', '-r300', '-nocrop')

figure(14)
grid on
plot(t, y_w(:,4))
xlabel('Time [s]')
ylabel('rb/2V_w [rad]')
set(gcf,'color','white')
legend('rb/2V_w')
export_fig('timedomain_rb_w', '-png', '-r300', '-nocrop')

figure(15)
grid on
plot(t, y_w(:,5))
xlabel('Time [s]')
ylabel('Lateral acceleration - a_{y_w} [m/s^2]')
set(gcf,'color','white')
legend('a_{y_w}')
export_fig('timedomain_ay_w', '-png', '-r300', '-nocrop')
