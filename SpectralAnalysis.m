% Filename : SpectralAnalysis (originally exampl83.m)
% From lecture notes ae4-304
%
% Bo Lee (5225604)
%%

% RUN PREVIOUS SCRIPTS
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
Pbeta  = (1/T)*( BETA.*conj(BETA));
Pphi   = (1/T)*(  PHI.*conj(PHI));
Pp     = (1/T)*(    P.*conj(P));
Pr     = (1/T)*(    R.*conj(R));
Pa_y   = (1/T)*(  A_Y.*conj(A_Y)); 

P_all = [Pbeta Pphi Pp Pr Pa_y];

Pbeta_r  = (1/T)*( BETA.*conj(BETA_r));
Pr_r     = (1/T)*(    R.*conj(R_r));
Pa_y_r   = (1/T)*(  A_Y.*conj(A_Y_r)); 

P_all_r = [Pbeta_r Pr_r Pa_y_r];

% DEFINE FREQUENCY VECTOR
fs = 1/dt;                                  % sample frequency
omega = 2*pi*fs*(0:(N/2)-1)/N;

% SMOOTHING FILTER
Pbeta_s  = zeros((length(Pbeta)),1);
Pphi_s   = zeros((length(Pphi)),1);
Pp_s     = zeros((length(Pp)),1);
Pr_s     = zeros((length(Pr)),1);
Pa_y_s  = zeros((length(Pa_y)),1);

P_all_s = [Pbeta_s Pphi_s Pp_s Pr_s Pa_y_s];

for ii = 1:length(P_all(1,:))
    for jj = 2:length(Pbeta)-2
        P_all_s(jj-1,ii) = 0.25*P_all(jj-1,ii)+0.5*P_all(jj,ii)+0.25*P_all(jj+1,ii);
    end
end

Pbeta_s_r  = zeros((length(Pbeta_r)),1);
Pr_s_r     = zeros((length(Pr_r)),1);
Pa_y_s_r  = zeros((length(Pa_y_r)),1);

P_all_s_r = [Pbeta_s_r Pr_s_r Pa_y_s_r];

for ii = 1:length(Periodograms_r(1,:))
    for jj = 2:length(Pbetar)-2
        P_all_s_r(jj-1,ii) = 0.25*P_all_r(jj-1,ii)+0.5*P_all_r(jj,ii)+0.25*P_all_r(jj+1,ii);
    end
end

% PLOTS
figure(14)
loglog(omega,Pbeta (1:N/2),'-',omega(1:length(omega)-1),P_all_s(1:N/2-1,1),w,Sxx(:,1),'k') 
axis(10.^[-1 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('Controlled S_\beta_\beta[rad^2/rad/s]');
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on
figure(15)
loglog(omega,Pphi(1:N/2),'-',omega(1:length(omega)-1),P_all_s(1:N/2-1,2),w,Sxx(:,2),'k');
axis(10.^[-1 2 -12 0]); xlabel('omega [rad/s]'); ylabel('Controlled S_\phi_\phi [rad^2/rad/s]')
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on
figure(16)
loglog(omega,Pp(1:N/2),'-',omega(1:length(omega)-1),P_all_s(1:N/2-1,3),w,Sxx(:,3),'k');
axis(10.^[-1 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('Controlled S_p_p [rad^2/rad/s]')
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on
figure(17)
loglog(omega,Pr(1:N/2),'-',omega(1:length(omega)-1),P_all_s(1:N/2-1,4),w,Sxx(:,4),'k');
axis(10.^[-1 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('Controlled S_r_r [rad^2/rad/s]')
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on
figure(18)
loglog(omega,Pa_y(1:N/2),'-',omega(1:length(omega)-1),P_all_s(1:N/2-1,5),w,Sxx(:,5),'k');
axis(10.^[-1 2 -10 5]); xlabel('omega [rad/s]'); ylabel('Controlled S_a_a [(m/s^2)^2/rad/s]')
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on

figure(19)
bode(A_cl,B,C(1,:),D(1,:),4,w)
hold on
bode(A_cl,B,C(1,:),D(1,:),5,w)
legend('lateral turbulence','vertical turbulence')

figure(20)
bode(A_cl,B,C(4,:),D(4,:),4,w)
hold on
bode(A_cl,B,C(4,:),D(4,:),5,w)
legend('lateral turbulence','vertical turbulence')

figure(21)
loglog(omega,Pbeta_r (1:N/2),'-',omega(1:length(omega)-1),P_all_s_r(1:N/2-1,1),w,Sxx(:,1),'k') 
axis(10.^[-1 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('Reduced S_\beta_\beta[rad^2/rad/s]');
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on
figure(24)
loglog(omega,Pr_r(1:N/2),'-',omega(1:length(omega)-1),P_all_s_r(1:N/2-1,2),w,Sxx(:,2),'k');
axis(10.^[-1 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('Reduced S_r_r [rad^2/rad/s]')
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on
figure(25)
loglog(omega,Pa_y_r(1:N/2),'-',omega(1:length(omega)-1),P_all_s_r(1:N/2-1,3),w,Sxx(:,3),'k');
axis(10.^[-1 2 -10 5]); xlabel('omega [rad/s]'); ylabel('Reduced S_a_a [(m/s^2)^2/rad/s]')
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on
