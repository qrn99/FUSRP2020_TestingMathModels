function dy = timeDelayODE(t,y,params)
%TIMEDELAYODE %This function will serve as the systems of ODE to be solved in the time
% delay model

%The inputs will be the timepoint t, y will be an array such that y(1)=S,
%y(2)=J, y(3)=L, y(4)=G and y(5)=I, and params will be an array with 17 parameters  

% params(1)=k_js;
% params(2)=k_gj;
% params(3)=k_jl;
% params(4)=k_gl;
% params(5)=k_xg;
% params(6)=k_xgi;
% params(7)=eta;
% params(8)=beta;
% params(9)=gamma;
% params(10)=f_gj;
% params(11)=k_xi;
% params(12)=k_lambda;
% params(13)=S_initial; 
% params(14)=J_initial; 
% params(15)=L_initial; 
% params(16)=I_initial; 
% params(17)=tau;

%initialize global variables
global G_initial;
global J;

%We assume a steady state and set the initial production of glucose,Gprod0,
%equal to the insulin dependent consumption of glucose, represented by 
%params(6)*G_initial*params(16), plus the glucose dependent consumption of 
%glucose, represented by params(5)*G_initial.
Gprod0=params(6)*G_initial*params(16)+params(5)*G_initial;
%We perform the following to calculate the hepatic production of glucose as
%detailed in the paper https://doi.org/10.3389/fbioe.2020.00195
steady_state_term=(params(12))/Gprod0;
denom=steady_state_term+(y(4)-G_initial);
hepatic_production=params(12)/denom;

%the function time_delay takes in inputs of current time point t, an array
%consisting of values of J for timepoints less than t, and our time delay
%constant tau
%the output is the glucose amount sent from J to L
function time_delay = delay(t,J,tau)
%for negative tau values let J send glucose to L with no time delay    
if tau<0
    time_delay=J(round(t)+1);
else
    %if t is less than tau then at timepoint t, no glucose from J is
    %sent to L
    if t<tau
        time_delay=0;
    %if t is greater than tau then at timepoint t, the glucose from J at 
    %timepoint t-tau is sent to L    
    else
        if round(round(t)-tau+1)>size(J,2)
            print('oh no!!!!')
        end
        time_delay=J(round(round(t)-tau+1));
    end
end
end
%We set up the systems of differential equations as detailed in the 
%paper https://doi.org/10.3389/fbioe.2020.00195
dy(1)=-params(1)*y(1);
dy(2)=params(1)*y(1)-params(2)*y(2)-params(3)*y(2);
J(round(t)+1)=y(2);
dy(3)=params(3)*delay(t,J,params(17))-params(4)*y(3);
dy(4)=-(params(5)+params(6)*y(5))*y(4)+hepatic_production+params(7)*(params(2)*y(2)+params(4)*y(3));
dy(5)=(params(11)*params(16))*((params(8)^(params(9))+1)/(params(8)^(params(9))*(G_initial/(y(4)+params(10)*(params(2)*y(2)+params(4)*y(3))))^params(9)+1)-y(5)/params(16));
dy=dy';
end

