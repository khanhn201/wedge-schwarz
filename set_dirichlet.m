function [U,V,P,T] = set_dirichlet(Uo,Vo,Po,To,Mu,Mv,Mt,X,Y);

N1 = size(X,1);
E = size(X,2);
Nelx = 4;

U = 0 - 0*X;  %% Desired field at inflow
% U(:,E-Nelx+1:E, N1) = 1;
V = 0 + 0*X;  %% Desired field at inflow
P = 0 + 0*X;  %% Desired field at inflow
T = 0 + 0*X;  %% Desired field at inflow
T(:,E-Nelx+1:E, N1) = 1;

U = Mu.*Uo + (1-Mu).*U;   %% Old value in interior, new value on inflow
V = Mv.*Vo + (1-Mv).*V;
P = Mv.*Po + (1-Mv).*P;
T = Mt.*To + (1-Mt).*T;
% for i = 1:E
%         reshape(U(:,i,:),N1, N1)
% end
% pause

% se_mesh(X,Y,U,'U bdry'); pause(1); pause



























