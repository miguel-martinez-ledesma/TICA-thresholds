% Levenberg-Marquardt OPTIONS
% 	for 4 Parameters 
% -------------------------------------------------------------------------
% Author: Miguel Martínez Ledesma (University of Chile)  
% Email: miguel.martinez@ing.uchile.cl
% Date: 19 November 2017
% -------------------------------------------------------------------------


%% Initialize Options Struct

LM_configuration={};

%% Algorithm Options
%             parameter    defaults    meaning
% opts(1)  =  prnt            3        >1 intermediate results; >2 plots; <0 no warning
% opts(2)  =  MaxIter      10*Npar     maximum number of iterations
% opts(3)  =  epsilon_1       1e-3     convergence tolerance for gradient
% opts(4)  =  epsilon_2       1e-3     convergence tolerance for parameters
% opts(5)  =  epsilon_3       1e-3     convergence tolerance for Chi-square
% opts(6)  =  epsilon_4       1e-1     determines acceptance of a L-M step
% opts(7)  =  lambda_0        1e-2     initial value of L-M paramter
% opts(8)  =  lambda_UP_fac   11       factor for increasing lambda
% opts(9)  =  lambda_DN_fac    9       factor for decreasing lambda
% opts(10) =  Update_Type      1       1: Levenberg-Marquardt lambda update
%                                      2: Quadratic update 
%                                      3: Nielsen's lambda update equations
% opts(11) =  MaxUnaccepted   500      Maximum number of consecutive unaccepted iterations
%
%                            prnt MaxIter   eps1    eps2   eps3   eps4 lam0  lamUP lamDN UpdateType MaxUnaccepted
LM_configuration.options  = [ -1,    200, 1e-12,  1e-18,     0, 1e-12, 1e-9,    11,    9,        1, 200]; 
% Note: lamda_0 changes the number of correct results and the probability to get locked in a local minimum for n>1e11

% LM weight 
% weighting vector for leasts quares fit ( weight >= 0 ) ... inverse of the standard measurement errors
% Default: sqrt (d.o.f. / (y dat’ *y dat))
LM_configuration.weight = 1;

%LM fractional increment of 'p' for numerical derivatives
LM_configuration.fractional_increment_p = -0.01;

%% Define Parameters

LM_configuration.num_var = 4;

%% Define Search Range

%                           [ Te    Ti    Vi    p ]
LM_configuration.par_min =  [ 200   200  -1000  0 ]; % minimum expected parameter values
LM_configuration.par_max =  [ 6000  6000   1000  1 ]; % maximum expected parameter values

