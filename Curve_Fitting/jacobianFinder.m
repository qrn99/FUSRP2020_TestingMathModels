function [jacob] = jacobianFinder(func,params)
%JACOBIANFINDER A numerical attempt to calculate the jacobian of one of our
%systems of ode; all implementations have failed thus far

%func is to represent our ode system, params are the variable parameters
%for which the jacobian is being calculated

    n=length(params);
    fevaluated=func(params);
    eps=1.e-8;
    paramperturb=params;
    jacob=zeros(length(fx),n);
    for i=1:n
        paramperturb(i)=paramperturb(i)+eps;
        jacob(:,i)=(func(paramperturb)-fevaluated)/eps;
        paramperturb(i)=params(i);
    end
end

