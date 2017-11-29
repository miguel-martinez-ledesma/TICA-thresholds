% y_hat = lm_func(t,p,c)
%  
% Function for nonlinear least squares curve-fitting
% using the Levenberg-Marquardt function, lm.m
%
%-----------------------------------------------------------------------------
% INPUT VARIABLES 
%  t     = m-vector of independent variable values (assumed to be error-free)
%  p     = n-vector of parameter values
%  c     = optional vector of other constants
%
%-----------------------------------------------------------------------------
% OUTPUT VARIABLES 
% y_hat  = m-vector of the curve-fit function evaluated at points t and
%          with parameters p
%  t     = m-vector of independent variable values (assumed to be error-free)
%
%-----------------------------------------------------------------------------
% Original from: 
%  H.P. Gavin, Dept. Civil & Environ. Eng'g, Duke Univ.
%
% Modified by: 
%	Miguel Martínez Ledesma (Universidad de Chile)
%	miguel.martinez@ing.uchile.cl
% Date: 19 November 2017
%-----------------------------------------------------------------------------



function [y_hat,t] = lm_func(t,p,c)

[y_hat,t] = run_ISR_simulator(c.radarFreq,c.NumPoints,c.maxFreq,c.Te,p(2,:),10^(p(1,:)),p(3,:),p(4,:));

% LM_FUNC ------------------------------------ 
