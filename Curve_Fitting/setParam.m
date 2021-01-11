function setParam
%Define global parameters with the following initial guesses  
global G_initial;
G_initial=5.9;%fasting glucose level mmol/L
global I_initial;
I_initial=1.7*10^(-7);%fasting insulin level mmol/L
global S_initial;
S_initial=417; %glucose ingested (mmol)
global J_initial; %jejunum
J_initial=0;
global L_initial; %ileum
L_initial=0;
global C_initial; %delay compartment
C_initial=0;
global k_js;
k_js=0.0664878906;%kinetic constant for glucose deliver from S to J
global k_gj;
k_gj=0.0788518734;%kinetic constant for glucose absorption from J
global k_jl;
k_jl=0.376187572;%kinetic constant for glucose deliver from J to L
global k_xg;
k_xg=0.01;%kinetic constant of basal uptake
global k_xgi;
k_xgi=10^(-6.5);%kinetic constant of insulin sensitive uptake
global k_gl;
k_gl=0.110271265; %kinetic constant for glucose absorption from L
global eta;
eta=0.01; %bioavaliability of absorbed glucose
global k_lambda
k_lambda=0.0316; %kinetic constant of glucose release
global f_gj
f_gj=10; %incretin action factor
global k_xi
k_xi=0.01; %kinetic constant for insulin degradatio
global beta
beta=100; %scale for insulin production saturation
global gamma
gamma=5; %scale for insulin production acceleration
global J
J=zeros(1,120);
global tau %time delay constant
tau=80;
global k_jc; %kinetic constant for gluocse delivery from J to C
k_jc=0.09;
global k_cl; %kinetic constant for gluocse delivery from C to L
k_cl=0.06;
end
