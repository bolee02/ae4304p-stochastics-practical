% Filename : Variances
% From lecture notes ae4-304
%
% Bo Lee (5225604)
%%

% RUN PREVIOUS SCRIPTS
SpectralAnalysis;
close all;

dw = diff(w);
dw(length(dw)+1) = 0;

var__beta_c = sum(Sxx(:,1)'.*dw)/pi;
var__phi_c = sum(Sxx(:,2)'.*dw)/pi;
var__pbv_c = sum(Sxx(:,3)'.*dw)/pi;
var__rbv_c = sum(Sxx(:,4)'.*dw)/pi;
var__ay_c = sum(Sxx(:,5)'.*dw)/pi;

ana = [var__beta_c var__phi_c var__pbv_c var__rbv_c var__ay_c];

vbeta_c  = var(beta_c);
vphi_c  = var(yt1(:,2));
vpbv_c  = var(yt1(:,3));
vrbv_c  = var(yt1(:,4));
vay_c  = var(a_y);

var_f = [vbeta_c vphi_c vpbv_c vrbv_c vay_c ];