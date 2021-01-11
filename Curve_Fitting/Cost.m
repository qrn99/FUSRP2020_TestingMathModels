function cost = Cost(params, gdata, version)
%COST: Depending on the model used, the cost function will output the sum of
%the square of the residuals
%   Version=1 corresponds to the no time delay model, 
%   Version=2 corresponds to the time delay model with time delay tau
%   Version=3 corresponds to the delay compartment model
%   gdata will be the collected glucose data we wish to fit our model to
%   params will be an array containing the value of our different
%   parameters for each glucose model; the number of parameters depends on
%   our chosen model, 16 for no time delay, 17 for time delay and 18 for
%   delay compartment

%   For each version, we will solve the systems of ODEs and take the
%   residual of predicted glucose values with the actual glucose data
%   Each residual wil be squared and the sum of the squared residuals will
%   be returned

%initialize global variables
global J;
global G_initial;

%No Time Delay Cost Function
if version==1
    %Our initial glucose value is given by the first data point of glucose 
    G_initial=gdata(1); 
    %Set the initial values of S,J,L,G,I
    init_val=[params(13) params(14) params(15) gdata(1) params(16)];
    %set the time span that we will solve our model for
    tspan = 0:(size(gdata)-1);
    %solve the chosen model's systems of ODE with the given initial conditions
    [~,y] = ode45(@(t,y) noTimeDelayODE(t,y,params),tspan,init_val);
    %cost is equal to the sum of the squares of the residuals
    cost=sum((gdata-y(:,4)).^2);
    %We add this penalty to ensure our minimization function does not
    %converge to a negative value for any of the parameters, since a
    %negative value is physiological impossible
    if sum(params>=0)~=16
        cost=cost+100;
    end
    
%Time Delay Cost Function
elseif version==2
    %Our initial glucose value is given by the first data point of glucose 
    G_initial=gdata(1);
    %The time delay component of this model requires us to store the value
    %of J at each calculated time point; we set variable J equal to the following
    %array of zeros so that as we calculate the value of J at different
    %time points, we store the values in variable J
    J=zeros(1,size(gdata,1)+2);
    %Set the initial values of S,J,L,G,I
    init_val=[params(13) params(14) params(15) gdata(1) params(16)];
    %set the time span that we will solve our model for
    tspan = 0:(size(gdata)-1);
    %solve the chosen model's systems of ODE with the given initial conditions
    [~,y] = ode45(@(t,y) timeDelayODE(t,y,params),tspan,init_val);
    %cost is equal to the sum of the squares of the residuals
    residuals=gdata-y(:,4);
    cost=sum(residuals.^2);
    % Penalty function as described earlier to ensure parameter values are
    % not negative
    if sum(params>=0)~=17
        cost=cost+100;
    end
    
%Delay Compartment Cost Function
elseif version==3
    %Our initial glucose value is given by the first data point of glucose 
    G_initial=gdata(1);
    %Set the initial values of S,J,L,G,I
    init_val=[params(12) params(13) params(14) gdata(1) params(15) params(18)];
    %set the time span that we will solve our model for
    tspan = 0:(size(gdata)-1);
    %solve the chosen model's systems of ODE with the given initial conditions
    [~,y] = ode45(@(t,y) delayCompartmentODE(t,y,params),tspan,init_val);
    residuals=gdata-y(:,4);
    %cost is equal to the sum of the squares of the residuals
    cost=sum(residuals.^2); 
    % Penalty function as described earlier to ensure parameter values are
    % not negative
    if sum(params>=0)~=18
        cost=cost+100;
    end
end
end

