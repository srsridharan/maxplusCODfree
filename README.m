%% README.m   
%% Purpose of the code/directory

% This directory contains the code for the reduced dimensionality approach to solving the HJB PDE for a quantum spin 
% system with unitary evolution

%% Dependency
% the cvx toolbox and the quantum optics toolbox have been used. Both of these are open source packages
% The latter can be worked around by simply recreating the Pauli matrices and the kron tensor.
% The cvx optimization routines are essential to the pruning although YALMIP or alternatives can also be used. 

%% Using the code
% The only information that needs to be entered is in 'MainFile.m' regarding the SMTP server, From and To email ids
% for notification of the results


%% Classes designed
% There are two classes : one for the codfreeobj (curse of dimensionality free object) class 
%		 This takes care of creating the codfree object, growing the control set  and using the CODfree approaches to
%		 pruning

%  The other class created is the plotcodfreeobj class.
%		 This is derived from the handle class and has routines to plot the performance of the curse of dimensionality 
%		 free approaches. It also has routines for error analysis.



