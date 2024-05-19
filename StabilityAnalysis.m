% Filename:  StabilityAnalysis.m (originally cita2a.m)
% From lecture notes ae4-304
%
% Bo Lee (5225604)
%%

% AIRCRAFT- AND FLIGHT CONDITION 'PRACTICAL'.
V   = 181.9;
S   = 24.2;
b   = 13.36;
mub = 32;
KX2 = 0.013;
KZ2 = 0.037;
KXZ = 0.002;
CL  = 1.1360;

% TURBULENCE PARAMETERS APPROXIMATED POWER SPECTRAL DENSITIES
Lg        = 150; 
B         = b/(2*Lg);
sigma     = 2;
sigmaug   = 3;
sigmaug_V = sigmaug/V;
sigmabg   = sigma/V;
sigmaag   = sigma/V;

Iug0 = 0.0249*sigmaug_V^2;
Iag0 = 0.0182*sigmaag^2;
tau1 = 0.0991;     tau2 = 0.5545;     tau3 = 0.4159;
tau4 = 0.0600;     tau5 = 0.3294;     tau6 = 0.2243;

% AIRCRAFT ASYMMETRIC AERODYNAMIC DERIVATIVES 
CYb  =-1.4900;     Clb  =-0.1240;     Cnb  = 0.1865;
CYp  =-0.1450;     Clp  =-0.4344;     Cnp  = 0.0135;
CYr  = 0.4300;     Clr  = 0.1550;     Cnr  =-0.1930;
CYda = 0.0000;     Clda =-0.2108;     Cnda = 0.0031;
CYdr = 0.3037;     Cldr = 0.0469;     Cndr =-0.1261;
 
                   Clpw = 0.8*Clp;    Cnpw = 0.9*Cnp;
                   Clrw = 0.7*Clr;    Cnrw = 0.2*Cnr;
CYfb = 0;
Clfb = 0;
Cnfb = 0;

%CYfbg = CYfb+0.5*CYr;
%Clfbg = Clfb+0.5*Clr;
%Cnfbg = Cnfb+0.5*Cnr;

% CALCULATION OF AIRCRAFT ASYMMETRIC STABILITY DERIVATIVES
yb   = (V/b)*CYb/(2*mub);
yphi = (V/b)*CL/(2*mub);
yp   = (V/b)*CYp/(2*mub);
yr   = (V/b)*(CYr-4*mub)/(2*mub);
ybg  = yb;
ydr  = (V/b)*CYdr/(2*mub);
den  = b*4*mub*(KX2*KZ2-KXZ^2)/V;
lb   = (Clb*KZ2+Cnb*KXZ)/den;
lp   = (Clp*KZ2+Cnp*KXZ)/den;
lr   = (Clr*KZ2+Cnr*KXZ)/den;
lda  = (Clda*KZ2+Cnda*KXZ)/den;
ldr  = (Cldr*KZ2+Cndr*KXZ)/den;
lug  = (-Clrw*KZ2-Cnrw*KXZ)/den;
lbg  = lb;
lag  = (Clpw*KZ2+Cnpw*KXZ)/den;
nb   = (Clb*KXZ+Cnb*KX2)/den;
np   = (Clp*KXZ+Cnp*KX2)/den;
nr   = (Clr*KXZ+Cnr*KX2)/den;
nda  = (Clda*KXZ+Cnda*KX2)/den;
ndr  = (Cldr*KXZ+Cndr*KX2)/den;
nug  = (-Clrw*KXZ-Cnrw*KX2)/den;
nbg  = nb;
nag  = (Clpw*KXZ+Cnpw*KX2)/den;
aug1 =-(V/Lg)^2*(1/(tau1*tau2));
aug2 =-(tau1+tau2)*(V/Lg)/(tau1*tau2);
aag1 =-(V/Lg)^2*(1/(tau4*tau5));
aag2 =-(tau4+tau5)*(V/Lg)/(tau4*tau5);
abg1 =-(V/Lg)^2;
abg2 =-2*(V/Lg);
bug1 = tau3*sqrt(Iug0*V/Lg)/(tau1*tau2);
bug2 = (1-tau3*(tau1+tau2)/(tau1*tau2))*sqrt(Iug0*(V/Lg)^3)/(tau1*tau2);
bag1 = tau6*sqrt(Iag0*V/Lg)/(tau4*tau5);
bag2 = (1-tau6*(tau4+tau5)/(tau4*tau5))*sqrt(Iag0*(V/Lg)^3)/(tau4*tau5);
bbg1 = sigmabg*sqrt(3*V/Lg);
bbg2 = (1-2*sqrt(3))*sigmabg*sqrt((V/Lg)^3);

% STATE- AND INPUT MATRICES
A = [yb yphi yp    yr 0    0    0    0    ybg  0;
     0  0    2*V/b 0  0    0    0    0    0    0;
     lb 0    lp    lr lug  0    lag  0    lbg  0;
     nb 0    np    nr nug  0    nag  0    nbg  0;
     0  0    0     0  0    1    0    0    0    0;
     0  0    0     0  aug1 aug2 0    0    0    0;
     0  0    0     0  0    0    0    1    0    0;
     0  0    0     0  0    0    aag1 aag2 0    0;
     0  0    0     0  0    0    0    0    0    1;
     0  0    0     0  0    0    0    0    abg1 abg2];

B = [0   ydr 0    0    0;
     0   0   0    0    0;
     lda ldr 0    0    0;
     nda ndr 0    0    0;
     0   0   bug1 0    0;
     0   0   bug2 0    0;
     0   0   0    bag1 0;
     0   0   0    bag2 0;
     0   0   0    0    bbg1;
     0   0   0    0    bbg2];

 C = eye(10,10);

 D = zeros(10,5);

% SHOW EIGENVALUES OF THE UNCONTROLLED SYSTEM
eig(A)
sys_ol = ss(A,B,C,D);
figure(1)
pzmap(sys_ol)
grid on

% SHOW EIGENVALUES OF THE REDUCED SYSTEM
A_r = [yb  yr 0    0    0    0    ybg  0;
       nb  nr nug  0    nag  0    nbg  0;
       0   0  0    1    0    0    0    0;
       0   0  aug1 aug2 0    0    0    0;
       0   0  0    0    0    1    0    0;
       0   0  0    0    aag1 aag2 0    0;
       0   0  0    0    0    0    0    1;
       0   0  0    0    0    0    abg1 abg2];

B_r = [0   ydr 0    0    0;
       nda ndr 0    0    0;
       0   0   bug1 0    0;
       0   0   bug2 0    0;
       0   0   0    bag1 0;
       0   0   0    bag2 0;
       0   0   0    0    bbg1;
       0   0   0    0    bbg2];

C_r = eye (8,8);

D_r = zeros(8,5); 

eig(A_r)
sys_r = ss(A_r, B_r, C_r, D_r);

figure(2)
pzmap(sys_r)
grid on

% SHOW EIGENVALUES OF THE CONTROLLED SYSTEM
Kphi = -0.1;
K    = [0 Kphi 0 0  0 0  0 0  0 0];
A_cl   = A-B(:,1)*K;
eig(A_cl)
sys_cl = ss(A_cl, B, C, D);

figure(3)
pzmap(sys_cl)
grid on