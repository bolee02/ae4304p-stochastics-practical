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

% ANALYTICAL POWER SPECTRA VARIANCE
var__beta_c = sum(Sxx(:,1)'.*dw)/pi;
var__phi_c = sum(Sxx(:,2)'.*dw)/pi;
var__pbv_c = sum(Sxx(:,3)'.*dw)/pi;
var__rbv_c = sum(Sxx(:,4)'.*dw)/pi;
var__ay_c = sum(Sxx(:,5)'.*dw)/pi;

ana = [var__beta_c var__phi_c var__pbv_c var__rbv_c var__ay_c];

var__beta_r = sum(Sxx_r(:,1)'.*dw)/pi;
var__rbv_r = sum(Sxx_r(:,2)'.*dw)/pi;
var__ay_r = sum(Sxx_r(:,3)'.*dw)/pi;

ana_r = [var__beta_r var__rbv_r var__ay_r];

% EXPERIMENTAL POWER SPECTRA VARIANCE
var__beta_c_e = sum(Pbeta'.*dw)/pi;
var__phi_c_e = sum(Pphi(:,2)'.*dw)/pi;
var__pbv_c_e = sum(Pp(:,3)'.*dw)/pi;
var__rbv_c_e = sum(Pr(:,4)'.*dw)/pi;
var__ay_c_e = sum(Pa_y(:,5)'.*dw)/pi;

var__beta_c_es = sum(Pbeta_s(:,1)'.*dw)/pi;
var__rbv_c_es = sum(Pr_s(:,2)'.*dw)/pi;
var__ay_c_es = sum(Pa_y_s(:,3)'.*dw)/pi;

% VARIANCE FROM var.m
vbeta_c  = var(beta_c);
vphi_c  = var(y(:,2));
vpbv_c  = var(y(:,3));
vrbv_c  = var(y(:,4));
vay_c  = var(a_y);
var_f = [vbeta_c vphi_c vpbv_c vrbv_c vay_c ];

vbeta_r  = var(ytr1(:,1));
vrbv_r  = var(ytr1(:,2));
vay_r  = var(a_y_r);
var_r = [vbeta_r vrbv_r vay_r];