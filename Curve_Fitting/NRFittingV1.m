setParam;
params = zeros(1,16);
global k_js
global k_gj
global k_jl
global k_xg
global k_xgi
global k_gl
global eta
global f_gj
global k_xi
global beta
global gamma
global k_lambda
global I_initial
global S_initial
global L_initial
global J_initial
params(1)=k_js;
params(2)=k_gj;
params(3)=k_jl;
params(4)=k_gl;
params(5)=k_xg;
params(6)=k_xgi;
params(7)=eta;
params(8)=beta;
params(9)=gamma;
params(10)=f_gj;
params(11)=k_xi;
params(12)=k_lambda;
params(13)=S_initial; %the initial glucose in stomach
params(14)=J_initial; %initial glucose in jejunum
params(15)=L_initial; %initiial glucose in ileum
params(16)=I_initial; %initial insulin concentration
tspan = 0:120;%120 time points, each representing a minute
init_val=[params(13) params(14) params(15) 5.995 params(16)];
data=readmatrix('2020-07-2021_lunch_peak.csv');
gdata=data(1:121,3);
fun=@(params)(Cost(params, gdata));
%x=lsqnonlin(fun,params)
x=fminsearch(fun,params)