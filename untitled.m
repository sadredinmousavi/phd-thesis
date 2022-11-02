function [x,fval,exitflag,output,lambda,grad,hessian] = untitled(Aineq,bineq,Aeq,beq,lb,ub)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimoptions('fmincon');
%% Modify options setting
options = optimoptions(options,'Display', 'off');
[x,fval,exitflag,output,lambda,grad,hessian] = ...
fmincon([],[],Aineq,bineq,Aeq,beq,lb,ub,@asfs,options);
