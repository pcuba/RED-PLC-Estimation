% Function to Construct GAUSS-HERMITE Integration Nodes
% This function returns the integration nodes. 
% 
% Updated (PC): 01/30/2017
%==========================================================================

function [EZ_nodes, WEIGHT] = get_GH_nodes(Q_in,par)


% quadrature nodes and weights
[quadpoint,weight] = ghpoints(Q_in);    

% Adjust quadrature points by stdev 
ea_prime = (sqrt(2)*par.STD_U*quadpoint);                    

% All the ez'
EZ_nodes = ea_prime;    

% weights for the product rule integral
WEIGHT  = ((1/sqrt(pi)))*(weight);  