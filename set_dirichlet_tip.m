function [U,V,P,T] = set_dirichlet_tip(Uo,Vo,Po,To,Mu,Mv,Mt,X,Y,omegax,omegay);

##U = 1 - Y.*Y;  %% Desired field at inflow
##V = 0 + 0*X;  %% Desired field at inflow
##T = 1 + 0*X;  %% Desired field at inflow

U = omegax;%0 + 0*X;  %% Desired field at inflow
V = omegay;%0 + 0*X;  %% Desired field at inflow
T = 0 + 0*X;  %% Desired field at inflow
P = 0 + 0*X;

U = Mu.*Uo + (1-Mu).*U;   %% Old value in interior, new value on inflow
V = Mv.*Vo + (1-Mv).*V;
P = Mv.*Po + (1-Mv).*P;
T = Mt.*To + (1-Mt).*T;

% se_mesh(X,Y,U,'U bdry'); pause(1); pause



























