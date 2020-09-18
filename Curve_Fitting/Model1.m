function Model1()
%MODEL1 Test the No Time Delay Model with 4 breakfast glucose curves and 4
%lunch glucose curves
%   The data on BreakfastToFit.xlsx is formatted such that each column
%   represents that data after eating breakfast on a particular day and
%   each row is a time point
%   Data from the BreakfastToFit excel sheets will read into the matrix data
%   such that the columns of data correspond to a particular day and the
%   rows to the appropriate time point. For each column, a curve fitting
%   algorithm will be employed to find the parameters that best allow our
%   model to predict the glucose data
%   The data on LunchToFit.xlsx is formatted such that each column
%   represents that data after eating lunch on a particular day and
%   each row is a time point
%   Data from the LunchToFit excel sheet will be read into the matrix data
%   just as before and the same proces will occur with this lunch data

setParam;
%initialize global variables
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

% test no time delay model
%set value of params
params = zeros(1,16);
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

%read Breakfast excel sheet into variable matrix data such that each column
%of data is glcuose data from a particular day and the rows are timepoints
data=readmatrix('BreakfastToFit.xlsx');

%gstore will be a variable that stores the actual glucose values in column 
%1 and the predicted glucose values from our model in column 2; we do not
%use gstore in our code but it is helpful for debugging
gstore=zeros(0,2);

%Notice size(data,2) reflects the number of days or glucose curves that we
%are fitting for

%storeBreakfastParam will be a 4x16 matrix storing the 16 parameter
%values in a designated column for each day as designated by the row
storeBreakfastParam=zeros(size(data,2),size(params,2));

%for each fitting we perform the following
for i=1:size(data,2)
    %we let gdata represent the actual glucose data we fit to
    gdata=rmmissing(data(:,i))*0.055; %rmmissing lets us ignore any missing row entries in the excel sheet after the final timepoint 
    %we set the maximum iterations of fminsearch to be 600 to prevent
    %excessive runtime
    options = optimset('MaxIter',600);
    %fun will now equal the current cost function, calculated as the sum of
    %the squared residuals
    fun=@(params)(Cost(params, gdata,1));
    %run the fminsearch algorithm minimizing the cost function value fun
    %while trying to optimize the parameters in params
    x= fminsearch(fun,params);
    %Set the initial values of S,J,L,G,I
    init_val=[x(13) x(14) x(15) gdata(1) x(16)];
    %set the time span that we will solve our model for
    tspan = 0:(size(gdata)-1);
    %solve the chosen model's systems of ODE with the given initial conditions
    [t,y] = ode45(@(t,y) noTimeDelayODE(t,y,x),tspan,init_val);
    
    %The following code comented out let's us graph the actual and
    %predicted glucose curves if desired
    
%     figure;
%     plot(tspan,gdata,'-o',tspan,y(:,4),'-o')
%     xlabel('Time (min)');
%     ylabel('Glucose (mM)');
%     title("Breakfast Day " + i)
%     legend('G Actual','G Modeled')

    %store the actual and predicted glucose values in gstore
    gstore=[gstore;[gdata, y(:,4)]];
    %we set the appropriate row of storeBreakfastParam equal to the
    %optimized params derived from fminsearch
    storeBreakfastParam(i,:)=x;
end
%Plot the following graph to determine the correlation between predicted
%and actual glucose values
ActualvsModel(gstore(:,1),gstore(:,2))
title("No Time Delay Model Breakfast")

%read Lunch excel sheet into variable matrix data such that each column
%of data is glcuose data from a particular day and the rows are timepoints
data=readmatrix('LunchToFit.xlsx');

%gstore will be a variable that stores the actual glucose values in column 
%1 and the predicted glucose values from our model in column 2; we do not
%use gstore in our code but it is helpful for debugging
gstore=zeros(0,2);

%Notice size(data,2) reflects the number of days or glucose curves that we
%are fitting for

%storeLunchParam will be a 4x16 matrix storing the 16 parameter
%values in a designated column for each day as designated by the row
storeLunchParam=zeros(size(data,2),size(params,2));

%for each fitting we perform the following
for i=1:size(data,2)
    %we let gdata represent the actual glucose data we fit to
    gdata=rmmissing(data(:,i))*0.055;%rmmissing lets us ignore any missing row entries in the excel sheet after the final timepoint
    %we set the maximum iterations of fminsearch to be 600 to prevent
    %excessive runtime
    options = optimset('MaxIter',600);
    %fun will now equal the current cost function, calculated as the sum of
    %the squared residuals
    fun=@(params)(Cost(params, gdata,1));
    %run the fminsearch algorithm minimizing the cost function value fun
    %while trying to optimize the parameters in params
    x= fminsearch(fun,params);
    %Set the initial values of S,J,L,G,I
    init_val=[x(13) x(14) x(15) gdata(1) x(16)];
    %set the time span that we will solve our model for
    tspan = 0:(size(gdata)-1);
    %solve the chosen model's systems of ODE with the given initial conditions
    [t,y] = ode45(@(t,y) noTimeDelayODE(t,y,x),tspan,init_val);

    %The following code comented out let's us graph the actual and
    %predicted glucose curves if desired
    
%     figure;
%     plot(tspan,gdata,'-o',tspan,y(:,4),'-o')
%     xlabel('Time (min)');
%     ylabel('Glucose (mM)');
%     title("Lunch Day " + i)
%     legend('G Actual','G Modeled')

    %store the actual and predicted glucose values in gstore
    gstore=[gstore;[gdata, y(:,4)]];
    %we set the appropriate row of storeLunchParam equal to the
    %optimized params derived from fminsearch
    storeLunchParam(i,:)=x;
end
%Plot the following graph to determine the correlation between predicted
%and actual glucose values
ActualvsModel(gstore(:,1),gstore(:,2))
title("No Time Delay Model Lunch")

%write optimized parameters for each curve fitting into the following excel
%file
writematrix([storeBreakfastParam; storeLunchParam],'NoTimeDelayParameters.xlsx')
end

