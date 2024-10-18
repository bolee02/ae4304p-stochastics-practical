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
temp = bode(A_cl,B,C(1,:),D(1,:),3,w);
Sbeta_u = temp.*temp;
temp = bode(A_cl,B,C(2,:),D(2,:),3,w);
Sphi_u = temp.*temp;
temp = bode(A_cl,B,C(3,:),D(3,:),3,w);
Spp_u = temp.*temp;
temp = bode(A_cl,B,C(4,:),D(4,:),3,w);
Srr_u = temp.*temp;
temp = bode(A_cl,B,C(5,:),D(5,:),3,w);
Sa_y_u = temp.*temp;

Sxx_u = [Sbeta_u Sphi_u Spp_u Srr_u Sa_y_u];

% RESPONSE TO VERTICAL TURBULENCE
temp = bode(A_cl,B,C(1,:),D(1,:),4,w); 
Sbeta_w = temp.*temp;
temp = bode(A_cl,B,C(2,:),D(2,:),4,w);
Sphi_w = temp.*temp;
temp = bode(A_cl,B,C(3,:),D(3,:),4,w);
Spp_w = temp.*temp;
temp = bode(A_cl,B,C(4,:),D(4,:),4,w);
Srr_w = temp.*temp;
temp = bode(A_cl,B,C(5,:),D(5,:),4,w);
Sa_y_w = temp.*temp;

Sxx_w = [Sbeta_w Sphi_w Spp_w Srr_w Sa_y_w];

% COMPUTE PSDS USING TIME DOMAIN DATA
u1 = [nn' nn' u_g' nn'  nn'];     
u2 = [nn' nn' nn'  w_g' nn'];

% COMPUTE SYSTEM RESPONSE
y_u     = lsim(A_cl,B,C,D,u1,t);
y_w     = lsim(A_cl,B,C,D,u2,t);

beta_u  = y_u(:,1);
phi_u   = y_u(:,2);
pbtV_u  = y_u(:,3);
rbtV_u  = y_u(:,4);
a_y_u   = y_u(:,5);

beta_w  = y_w(:,1);
phi_w   = y_w(:,2);
pbtV_w  = y_w(:,3);
rbtV_w  = y_w(:,4);
a_y_w   = y_w(:,5);

% COMPUTE PERIODOGRAM AND ESTIMATE PSD
% PERIODOGRAM
BETA_u  = dt*fft(beta_u);
PHI_u   = dt*fft(phi_u);
P_u     = dt*fft(pbtV_u);
R_u     = dt*fft(rbtV_u);
A_Y_u   = dt*fft(a_y_u);

BETA_w  = dt*fft(beta_w);
PHI_w   = dt*fft(phi_w);
P_w     = dt*fft(pbtV_w);
R_w     = dt*fft(rbtV_w);
A_Y_w   = dt*fft(a_y_w);

% PSD ESTIMATE
N_half = floor(N/2);

Pbeta_u = (2/T) * abs(BETA_u(1:N_half)).^2;
Pphi_u = (2/T) * abs(PHI_u(1:N_half)).^2;
Pp_u = (2/T) * abs(P_u(1:N_half)).^2;
Pr_u = (2/T) * abs(R_u(1:N_half)).^2;
Pa_y_u = (2/T) * abs(A_Y_u(1:N_half)).^2;

P_all_u = [Pbeta_u Pphi_u Pp_u Pr_u Pa_y_u];

Pbeta_w = (2/T) * abs(BETA_w(1:N_half)).^2;
Pphi_w = (2/T) * abs(PHI_w(1:N_half)).^2;
Pp_w = (2/T) * abs(P_w(1:N_half)).^2;
Pr_w = (2/T) * abs(R_w(1:N_half)).^2;
Pa_y_w = (2/T) * abs(A_Y_w(1:N_half)).^2;

P_all_w = [Pbeta_w Pphi_w Pp_w Pr_w Pa_y_w];

% DEFINE FREQUENCY VECTOR
fs = 1/dt;                                  % sample frequency
omega = 2*pi*fs*(0:N_half-1)/N;

% SMOOTHING FILTER
filter = [0.25, 0.5, 0.25];

P_all_s_u = conv2(P_all_u, filter(:), 'valid');
P_all_s_w = conv2(P_all_w, filter(:), 'valid');

% PLOTS for u
figure(16)
grid on
loglog(omega, Pbeta_u(1:N_half), omega(1:length(omega)-2), P_all_s_u(1:N_half-2,1), w, Sxx_u(:,1), 'k--') 
axis(10.^[-1 2 -12 -2]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{\beta\beta} [rad^2/Hz]');
legend('Periodogram', 'Smoothed Periodogram', 'Analytical PSD')
set(gcf, 'color', 'white')
export_fig('spec_beta_u', '-png', '-r300', '-nocrop')

figure(17)
grid on
loglog(omega, Pphi_u(1:N_half), omega(1:length(omega)-2), P_all_s_u(1:N_half-2,2), w, Sxx_u(:,2), 'k--');
axis(10.^[-1 2 -10 0]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{\phi\phi} [rad^2/Hz]')
legend('Periodogram', 'Smoothed Periodogram', 'Analytical PSD')
set(gcf, 'color', 'white')
export_fig('spec_phi_u', '-png', '-r300', '-nocrop')

figure(18)
grid on
loglog(omega, Pp_u(1:N_half), omega(1:length(omega)-2), P_all_s_u(1:N_half-2,3), w, Sxx_u(:,3), 'k--');
axis(10.^[-1 2 -12 -2]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{pp} [rad^2/Hz]')
legend('Periodogram', 'Smoothed Periodogram', 'Analytical PSD')
set(gcf, 'color', 'white')
export_fig('spec_p_u', '-png', '-r300', '-nocrop')

figure(19)
grid on
loglog(omega, Pr_u(1:N_half), omega(1:length(omega)-2), P_all_s_u(1:N_half-2,4), w, Sxx_u(:,4), 'k--');
axis(10.^[-1 2 -14 -2]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{rr} [rad^2/Hz]')
legend('Periodogram', 'Smoothed Periodogram', 'Analytical PSD')
set(gcf, 'color', 'white')
export_fig('spec_r_u', '-png', '-r300', '-nocrop')

figure(20)
grid on
loglog(omega, Pa_y_u(1:N_half), omega(1:length(omega)-2), P_all_s_u(1:N_half-2,5), w, Sxx_u(:,5), 'k--');
axis(10.^[-1 2 -6 4]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{aa} [(m/s^2)^2/Hz]')
legend('Periodogram', 'Smoothed', 'Analytical')
set(gcf, 'color', 'white')
export_fig('spec_a_u', '-png', '-r300', '-nocrop')

% PLOTS for w
figure(23)
grid on
loglog(omega, Pbeta_w(1:N_half), omega(1:length(omega)-2), P_all_s_w(1:N_half-2,1), w, Sxx_w(:,1), 'k--')
axis(10.^[-1 2 -12 -2]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{\beta_w\beta_w} [rad^2/Hz]');
legend('Periodogram', 'Smoothed Periodogram', 'Analytical PSD')
set(gcf, 'color', 'white')
export_fig('spec_beta_w', '-png', '-r300', '-nocrop')

figure(24)
grid on
loglog(omega, Pphi_w(1:N_half), omega(1:length(omega)-2), P_all_s_w(1:N_half-2,2), w, Sxx_w(:,2), 'k--');
axis(10.^[-1 2 -10 0]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{\phi_w\phi_w} [rad^2/Hz]')
legend('Periodogram', 'Smoothed Periodogram', 'Analytical PSD')
set(gcf, 'color', 'white')
export_fig('spec_phi_w', '-png', '-r300', '-nocrop')

figure(25)
grid on
loglog(omega, Pp_w(1:N_half), omega(1:length(omega)-2), P_all_s_w(1:N_half-2,3), w, Sxx_w(:,3), 'k--');
axis(10.^[-1 2 -12 -2]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{pp} [rad^2/Hz]')
legend('Periodogram', 'Smoothed Periodogram', 'Analytical PSD')
set(gcf, 'color', 'white')
export_fig('spec_p_w', '-png', '-r300', '-nocrop')

figure(26)
grid on
loglog(omega, Pr_w(1:N_half), omega(1:length(omega)-2), P_all_s_w(1:N_half-2,4), w, Sxx_w(:,4), 'k--');
axis(10.^[-1 2 -14 -2]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{rr} [rad^2/Hz]')
legend('Periodogram', 'Smoothed Periodogram', 'Analytical PSD')
set(gcf, 'color', 'white')
export_fig('spec_r_w', '-png', '-r300', '-nocrop')

figure(27)
grid on
loglog(omega, Pa_y_w(1:N_half), omega(1:length(omega)-2), P_all_s_w(1:N_half-2,5), w, Sxx_w(:,5), 'k--');
axis(10.^[-1 2 -6 4]); xlabel('Frequency - \omega [rad/s]'); ylabel('Controlled S_{aa} [(m/s^2)^2/Hz]')
legend('Periodogram', 'Smoothed', 'Analytical PSD')
set(gcf, 'color', 'white')
export_fig('spec_a_w', '-png', '-r300', '-nocrop')

% BODE PLOTS
figure(28)
bode(A_cl,B,C(1,:),D(1,:),3,w)
hold on
bode(A_cl,B,C(1,:),D(1,:),4,w)
legend('Longitudinal Acceleration - u_g', 'Vertical Acceleration - w_g')
set(gcf, 'color', 'white')
export_fig('bode_1', '-png', '-r300', '-nocrop')

figure(29)
bode(A_cl,B,C(2,:),D(2,:),3,w)
hold on
bode(A_cl,B,C(2,:),D(2,:),4,w)
legend('Longitudinal Acceleration - u_g', 'Vertical Acceleration - w_g')
set(gcf, 'color', 'white')
export_fig('bode_2', '-png', '-r300', '-nocrop')

figure(30)
bode(A_cl,B,C(3,:),D(3,:),3,w)
hold on
bode(A_cl,B,C(3,:),D(3,:),4,w)
legend('Longitudinal Acceleration - u_g', 'Vertical Acceleration - w_g')
set(gcf, 'color', 'white')
export_fig('bode_3', '-png', '-r300', '-nocrop')

figure(31)
bode(A_cl,B,C(4,:),D(4,:),3,w)
hold on
bode(A_cl,B,C(4,:),D(4,:),4,w)
legend('Longitudinal Acceleration - u_g', 'Vertical Acceleration - w_g')
set(gcf, 'color', 'white')
export_fig('bode_4', '-png', '-r300', '-nocrop')

figure(32)
bode(A_cl,B,C(5,:),D(5,:),3,w)
hold on
bode(A_cl,B,C(5,:),D(5,:),4,w)
legend('Longitudinal Acceleration - u_g', 'Vertical Acceleration - w_g')
set(gcf, 'color', 'white')
export_fig('bode_5', '-png', '-r300', '-nocrop')