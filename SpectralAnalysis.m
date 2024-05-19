% Filename : SpectralAnalysis (originally exampl83.m)
% From lecture notes ae4-304
%
% Bo Lee (5225604)
%%

% RUN PREVIOUS SCRIPTS
clc, clf, clear
TimeDomainAnalysis;
close all;

% DEFINE FREQUENCY VECTOR
w = logspace(-2,2,300);

% COMPUTE ANALYTIC POWER SPECTRAL DENSITIES
% RESPONSE TO HORIZONTAL LATERAL TURBULENCE
temp = bode(A_cl,B,C(1,:),D(1,:),4,w)+bode(A_cl,B,C(1,:),D(1,:),5,w); Sbeta  = temp.*temp;
temp = bode(A_cl,B,C(2,:),D(2,:),4,w)+bode(A_cl,B,C(2,:),D(2,:),5,w); Sphi   = temp.*temp;
temp = bode(A_cl,B,C(3,:),D(3,:),4,w)+bode(A_cl,B,C(3,:),D(3,:),5,w); Spp    = temp.*temp;
temp = bode(A_cl,B,C(4,:),D(4,:),4,w)+bode(A_cl,B,C(4,:),D(4,:),5,w); Srr    = temp.*temp;
temp = bode(A_cl,B,C(5,:),D(5,:),4,w)+bode(A_cl,B,C(5,:),D(5,:),5,w); Sa_y   = temp.*temp;

Sxx  = [Sbeta Sphi Spp Srr Sa_y];

temp = bode(A_r,B_r,C_r(1,:),D_r(1,:),4,w)+bode(A_r,B_r,C_r(1,:),D_r(1,:),5,w); Sbeta_r  = temp.*temp;
temp = bode(A_r,B_r,C_r(2,:),D_r(2,:),4,w)+bode(A_r,B_r,C_r(2,:),D_r(2,:),5,w); Srr_r   = temp.*temp;
temp = bode(A_r,B_r,C_r(3,:),D_r(3,:),4,w)+bode(A_r,B_r,C_r(3,:),D_r(3,:),5,w); Sa_y_r    = temp.*temp;

Sxx_r = [Sbeta_r Srr_r Sa_y_r];

% COMPUTE PSDS USING TIME DOMAIN DATA
u = [nn' nn' u_g' w_g'  nn'];     

% COMPUTE SYSTEM RESPONSE
y     = lsim(A_cl,B,C,D,u,t);
y_r     = lsim(A_r,B_r,C_r,D_r,u,t);

beta  = y(:,1);
phi   = y(:,2);
pbtV  = y(:,3);
rbtV  = y(:,4);
a_y   = y(:,5);

beta_r = y_r(:,1);
rbt_r  = y_r(:,2);
a_y_r  = y_r(:,3);

% COMPUTE PERIODOGRAM AND ESTIMATE PSD
% PERIODOGRAM
BETA  = dt*fft(beta);
PHI   = dt*fft(phi);
P     = dt*fft(pbtV);
R     = dt*fft(rbtV);
A_Y   = dt*fft(a_y);

BETA_r  = dt*fft(beta_r);
R_r     = dt*fft(rbtV);
A_Y_r   = dt*fft(a_y_r);

% PSD ESTIMATE
N_half = floor(N / 2);

Pbeta = (2/T) * abs(BETA(1:N_half)).^2 * 2*V/b;
Pphi = (2/T) * abs(PHI(1:N_half)).^2;
Pp = (2/T) * abs(P(1:N_half)).^2 * 2*V/b;
Pr = (2/T) * abs(R(1:N_half)).^2 * 2*V/b;
Pa_y = (2/T) * abs(A_Y(1:N_half)).^2 * 2*V/b;

P_all = [Pbeta Pphi Pp Pr Pa_y];

Pbeta_r = (2/T) * abs(BETA_r(1:N_half)).^2;
Pr_r = (2/T) * abs(R_r(1:N_half)).^2;
Pa_y_r = (2/T) * abs(A_Y_r(1:N_half)).^2;

P_all_r = [Pbeta_r Pr_r Pa_y_r];

% DEFINE FREQUENCY VECTOR
fs = 1/dt;                                  % sample frequency
omega = 2*pi*fs*(0:N_half-1)/N;

% SMOOTHING FILTER
filter = [0.25, 0.5, 0.25];

P_all_s = conv2(P_all, filter(:), 'valid');
P_all_s_r = conv2(P_all_r, filter(:), 'valid');

% PLOTS
figure(14)
grid on
loglog(omega,Pbeta(1:N/2),omega(1:length(omega)-2),P_all_s(1:N/2-2,1),w,Sxx(:,1),'k--') 
axis(10.^[-1 2 -12 -2]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{\beta\beta} [rad^2/rad/s]');
legend('Periodogram','Smoothed Periodogram','Analytical PSD')

figure(15)
grid on
loglog(omega,Pphi(1:N/2),omega(1:length(omega)-2),P_all_s(1:N/2-2,2),w,Sxx(:,2),'k--');
axis(10.^[-1 2 -10 0]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{\phi\phi} [rad^2/rad/s]')
legend('Periodogram','Smoothed Periodogram','Analytical PSD')

figure(16)
grid on
loglog(omega,Pp(1:N/2),omega(1:length(omega)-2),P_all_s(1:N/2-2,3),w,Sxx(:,3),'k--');
axis(10.^[-1 2 -12 -2]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{pp} [rad^2/rad/s]')
legend('Periodogram','Smoothed Periodogram','Analytical PSD')

figure(17)
grid on
loglog(omega,Pr(1:N/2),omega(1:length(omega)-2),P_all_s(1:N/2-2,4),w,Sxx(:,4),'k--');
axis(10.^[-1 2 -14 -2]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{rr} [rad^2/rad/s]')
legend('Periodogram','Smoothed Periodogram','Analytical PSD')

figure(18)
grid on
loglog(omega,Pa_y(1:N/2),omega(1:length(omega)-2),P_all_s(1:N/2-2,5),w,Sxx(:,5),'k--');
axis(10.^[-1 2 -6 4]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{aa} [(m/s^2)^2/rad/s]')
legend('Periodogram','Smoothed','Analytical')

figure(19)
grid on
loglog(omega,Pbeta_r(1:N/2),omega(1:length(omega)-2),P_all_s_r(1:N/2-2,1),w,Sxx(:,1),'k--') 
axis(10.^[-1 2 -12 -2]); xlabel('Frequency - \omega [rad/s]'); ylabel('Reduced S_{\beta\beta} [rad^2/rad/s]');
legend('Periodogram','Smoothed','Analytical')

figure(20)
grid on
loglog(omega,Pr_r(1:N/2),omega(1:length(omega)-2),P_all_s_r(1:N/2-2,2),w,Sxx(:,2),'k--');
axis(10.^[-1 2 -14 0]); xlabel('Frequency - \omega [rad/s]'); ylabel('Reduced S_{rr} [rad^2/rad/s]')
legend('Periodogram','Smoothed','Analytical')

figure(21)
grid on
loglog(omega,Pa_y_r(1:N/2),omega(1:length(omega)-2),P_all_s_r(1:N/2-2,3),w,Sxx(:,3),'k--');
axis(10.^[-1 2 -10 0]); xlabel('Frequency - \omega [rad/s]'); ylabel('Reduced S_{aa} [(m/s^2)^2/rad/s]')
legend('Periodogram','Smoothed','Analytical')

figure(22)
bode(A_cl,B,C(1,:),D(1,:),4,w)
hold on
bode(A_cl,B,C(1,:),D(1,:),5,w)
legend('Longitudinal Acceleration - u_g','Vertical Acceleration - w_g')

figure(23)
bode(A_cl,B,C(4,:),D(4,:),4,w)
hold on
bode(A_cl,B,C(4,:),D(4,:),5,w)
legend('Longitudinal Acceleration - u_g','Vertical Acceleration - w_g')