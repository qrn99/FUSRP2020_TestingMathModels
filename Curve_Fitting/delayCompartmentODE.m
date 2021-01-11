function dy = delayCompartmentODE(t,y,params)
%DELAYCOMPARTMENTODE % This function will serve as the systems of ODE to be
%solved in the compartment delay model

%The inputs will be the timepoint t, y will be an array such that y(1)=S, 
%y(2)=J, y(3)=L, y(4)=G, Y(5)=I, y(6)=C with C as the delayed compartment, 
%params will be an array with 18 parameters  

% params(1)=k_js;
% params(2)=k_gj;
% params(3)=k_gl;
% params(4)=k_xg;
% params(5)=k_xgi;
% params(6)=eta;
% params(7)=beta;
% params(8)=gamma;
% params(9)=f_gj;
% params(10)=k_xi;
% params(11)=k_lambda;
% params(12)=S_initial; 
% params(13)=J_initial; 
% params(14)=L_initial; 
% params(15)=I_initial; 
% params(16)=k_jc
% params(17)=k_cl
% params(18)=C_initial;

%initialize global variables
global G_initial;
%We assume a steady state and set the initial production of glucose,Gprod0,
%equal to the insulin dependent consumption of glucose, represented by 
%params(6)*G_initial*params(16), plus the glucose dependent consumption of 
%glucose, represented by params(5)*G_initial.
Gprod0=params(5)*G_initial*params(15)+params(4)*G_initial;
%We perform the following to calculate the hepatic production of glucose as
%detailed in the paper https://doi.org/10.3389/fbioe.2020.00195
steady_state_term=(params(11))/Gprod0;
denom=steady_state_term+(y(4)-G_initial);
hepatic_production=params(11)/denom;

%We set up the systems of differential equations as detailed in the 
%paper https://doi.org/10.3389/fbioe.2020.00195
%However, we set the time delay tau to 0
%Let J empty into C and C empty into L with kjc the rate costant of the
%former and kcl the rate constant of the latter; this informs dy(6)
dy(1)=-params(1)*y(1);
dy(2)=params(1)*y(1)-params(2)*y(2)-params(16)*y(2);
dy(3)=params(17)*y(6)-params(3)*y(3);
dy(4)=-(params(4)+params(5)*y(5))*y(4)+hepatic_production+params(6)*(params(2)*y(2)+params(3)*y(3));
dy(5)=(params(10)*params(15))*((params(7)^(params(8))+1)/(params(7)^(params(8))*(G_initial/(y(4)+params(9)*(params(2)*y(2)+params(3)*y(3))))^params(8)+1)-y(5)/params(15));
dy(6)= params(16)*y(2)-params(17)*y(6);
dy=dy';
end

