clear all;
setParam;
solveSystem;
function setParam
%Define parameters
global G_initial;
G_initial=5.8;%fasting glucose level mmol/L
global I_initial;
I_initial=1.7*10^(-7);%fasting insulin level units unknown and this is more of a guess mmol/L
global S_initial;
S_initial=417; %glucose ingested (mmol)
global J_initial;
J_initial=0;
global J;
J=zeros([1 121]);%this will store glucose at J for each minute; needed to time delay function
global L_initial;
L_initial=0;
global k_js;
k_js=0.1;%kinetic constant for glucose deliver from S to J
global k_gj;
k_gj=0.1;%kinetic constant for glucose absorption from J
global k_jl;
k_jl=0.0316;%kinetic constant for glucose deliver from J to L
global tau;
tau=80;%time delay constant between J to L
global k_xg;
k_xg=0.01;%kinetic constant of basal uptake
global k_xgi;
k_xgi=10^(-6.5);%kinetic constant of insulin sensitive uptake
global k_gl;
k_gl=0.1; %kinetic constant for glucose absorption from L
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
global Gprod0
Gprod0=k_xgi*G_initial*I_initial+k_xg*G_initial; %fasting_hepatic_contribution???????
end
function time_delay = delay(t,J)
global tau;
if t<tau
    time_delay=0;
else
    time_delay=J(round(round(t)-tau+1));
end
end
function hepatic_production = Gprod(G)
global k_lambda;
global Gprod0;
global G_initial;
    steady_state_term=(k_lambda)/Gprod0;
    denom=steady_state_term+(G-G_initial);
    hepatic_production=k_lambda/denom;
end

%let y(1)=S,y(2)=J, y(3)=L, y(4)=G, y(5)=I
function dy = myODE(t,y)
global G_initial
global I_initial
global k_js
global k_gj
global k_jl
global tau
global k_xg
global k_xgi
global k_gl
global eta
global f_gj
global k_xi
global beta
global gamma
global J;
dy(1)=-k_js*y(1);
dy(2)=k_js*y(1)-k_gj*y(2)-k_jl*y(2);
J(round(t)+1)=y(2);
dy(3)=k_jl*delay(t,J)-k_gl*y(3);
dy(4)=-(k_xg+k_xgi*y(5))*y(4)+Gprod(y(4))+eta*(k_gj*y(2)+k_gl*y(3));
dy(5)=(k_xi*I_initial)*((beta^(gamma)+1)/(beta^(gamma)*(G_initial/(y(4)+f_gj*(k_gj*y(2)+k_gl*y(3))))^gamma+1)-y(5)/I_initial);
dy=dy';
end
function solveSystem
global G_initial
global I_initial
global S_initial
global L_initial
global J_initial
tspan = 0:120;%120 time points, each representing a minute
init_val=[S_initial J_initial L_initial G_initial I_initial];
[t,y] = ode45(@myODE,tspan,init_val);
plot(t,y(:,1),'-o',t,y(:,2),'-o',t,y(:,3),'-o')
title('Solution with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('y(1)=S','y(2)=J', 'y(3)=L')
figure;
plot(t,y(:,4),'-o')
title('Solution with ODE45');
xlabel('Time t');
ylabel('Solution G');
legend('y(4)=G')
figure;
plot(t,y(:,5),'-o')
title('Solution with ODE45');
xlabel('Time t');
ylabel('Solution I');
legend('y(5)=I')
end
