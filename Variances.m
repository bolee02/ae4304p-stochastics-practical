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
var_beta_u = sum(Sxx_u(:,1)'.*dw)/pi;
var_phi_u = sum(Sxx_u(:,2)'.*dw)/pi;
var_pp_u = sum(Sxx_u(:,3)'.*dw)/pi;
var_rr_u = sum(Sxx_u(:,4)'.*dw)/pi;
var_a_y_u = sum(Sxx_u(:,5)'.*dw)/pi;

var_beta_w = sum(Sxx_w(:,1)'.*dw)/pi;
var_phi_w = sum(Sxx_w(:,2)'.*dw)/pi;
var_pp_w = sum(Sxx_w(:,3)'.*dw)/pi;
var_rr_w = sum(Sxx_w(:,4)'.*dw)/pi;
var_a_y_w = sum(Sxx_w(:,5)'.*dw)/pi;

% EXPERIMENTAL POWER SPECTRA VARIANCE
var_beta_c_e_u = sum(Pbeta_u(1:length(domega))'.*domega)/pi;
var_phi_c_e_u = sum(Pphi_u(1:length(domega))'.*domega)/pi;
var_pp_c_e_u = sum(Pp_u(1:length(domega))'.*domega)/pi;
var_rr_c_e_u = sum(Pr_u(1:length(domega))'.*domega)/pi;
var_a_y_c_e_u = sum(Pa_y_u(1:length(domega))'.*domega)/pi;

var_beta_c_e_w = sum(Pbeta_w(1:length(domega))'.*domega)/pi;
var_phi_c_e_w = sum(Pphi_w(1:length(domega))'.*domega)/pi;
var_pp_c_e_w = sum(Pp_w(1:length(domega))'.*domega)/pi;
var_rr_c_e_w = sum(Pr_w(1:length(domega))'.*domega)/pi;
var_a_y_c_e_w = sum(Pa_y_w(1:length(domega))'.*domega)/pi;

var_beta_c_es_u = sum(P_all_s_u(1:length(domega)-2,1)'.*domega(1:length(domega)-2))/pi;
var_phi_c_es_u = sum(P_all_s_u(1:length(domega)-2,2)'.*domega(1:length(domega)-2))/pi;
var_pp_c_es_u = sum(P_all_s_u(1:length(domega)-2,3)'.*domega(1:length(domega)-2))/pi;
var_rr_c_es_u = sum(P_all_s_u(1:length(domega)-2,4)'.*domega(1:length(domega)-2))/pi;
var_a_y_c_es_u = sum(P_all_s_u(1:length(domega)-2,5)'.*domega(1:length(domega)-2))/pi;

var_beta_c_es_w = sum(P_all_s_w(1:length(domega)-2,1)'.*domega(1:length(domega)-2))/pi;
var_phi_c_es_w = sum(P_all_s_w(1:length(domega)-2,2)'.*domega(1:length(domega)-2))/pi;
var_pp_c_es_w = sum(P_all_s_w(1:length(domega)-2,3)'.*domega(1:length(domega)-2))/pi;
var_rr_c_es_w = sum(P_all_s_w(1:length(domega)-2,4)'.*domega(1:length(domega)-2))/pi;
var_a_y_c_es_w = sum(P_all_s_w(1:length(domega)-2,5)'.*domega(1:length(domega)-2))/pi;

% VARIANCE FROM var.m
var_beta_f_u  = var(y_u(:,1));
var_phi_f_u  = var(y_u(:,2));
var_pbv_f_u  = var(y_u(:,3));
var_rbv_f_u  = var(y_u(:,4));
var_ay_f_u  = var(y_u(:,5));
var_f_u = [var_beta_f_u var_phi_f_u var_pbv_f_u var_rbv_f_u var_ay_f_u];

var_beta_f_w  = var(y_w(:,1));
var_phi_f_w  = var(y_w(:,2));
var_pbv_f_w  = var(y_w(:,3));
var_rbv_f_w  = var(y_w(:,4));
var_ay_f_w  = var(y_w(:,5));
var_f_w = [var_beta_f_w var_phi_f_w var_pbv_f_w var_rbv_f_w var_ay_f_w];
