% Filename : Variances
% From lecture notes ae4-304
%
% Bo Lee (5225604)
%%

% RUN PREVIOUS SCRIPTS
SpectralAnalysis;
close all;

dw = diff(w);
domega = diff(omega);

dw(length(dw)+1) = 0;
domega(length(domega)+1) = 0;

% ANALYTICAL POWER SPECTRA VARIANCE
var_beta = sum(Sxx(:,1)'.*dw)/pi;
var_phi = sum(Sxx(:,2)'.*dw)/pi;
var_pp = sum(Sxx(:,3)'.*dw)/pi;
var_rr = sum(Sxx(:,4)'.*dw)/pi;
var_a_y = sum(Sxx(:,5)'.*dw)/pi;

var_beta_r = sum(Sxx_r(:,1)'.*dw)/pi;
var_rr_r = sum(Sxx_r(:,2)'.*dw)/pi;
var_a_y_r = sum(Sxx_r(:,3)'.*dw)/pi;

% EXPERIMENTAL POWER SPECTRA VARIANCE
var_beta_c_e = sum(Pbeta(1:length(domega))'.*domega)/pi;
var_phi_c_e = sum(Pphi(1:length(domega))'.*domega)/pi;
var_pp_c_e = sum(Pp(1:length(domega))'.*domega)/pi;
var_rr_c_e = sum(Pr(1:length(domega))'.*domega)/pi;
var_a_y_c_e = sum(Pa_y(1:length(domega))'.*domega)/pi;

var_beta_c_es = sum(P_all_s(1:length(domega)-2,1)'.*domega(1:length(domega)-2))/pi;
var_phi_c_es = sum(P_all_s(1:length(domega)-2,2)'.*domega(1:length(domega)-2))/pi;
var_pp_c_es = sum(P_all_s(1:length(domega)-2,3)'.*domega(1:length(domega)-2))/pi;
var_rr_c_es = sum(P_all_s(1:length(domega)-2,4)'.*domega(1:length(domega)-2))/pi;
var_a_y_c_es = sum(P_all_s(1:length(domega)-2,5)'.*domega(1:length(domega)-2))/pi;

% VARIANCE FROM var.m
var_beta_f  = var(y(:,1));
var_phi_f  = var(y(:,2));
vpbv_c  = var(y(:,3));
vrbv_c  = var(y(:,4));
vay_c  = var(y(:,5));
var_f = [vbeta_c vphi_c vpbv_c vrbv_c vay_c ];

vbeta_r  = var(y_r(:,1));
vrbv_r  = var(y_r(:,2));
vay_r  = var(y_r(:,3));
var_r = [vbeta_r vrbv_r vay_r];