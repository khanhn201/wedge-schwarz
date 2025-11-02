function [JM,DM,BMh,Jf,dt] = set_dealiasing(N,z,Jac,U,V,Rx)
%% Set dealiasing operators, JM,DM,BMh
Nd = floor(1.5*N);
[zd,wd]=zwgl(Nd); Bd=diag(wd);
JM=interp_mat(zd,z);
DM=deriv_mat (zd);
BMh=tensor3(Bd*JM,1,Bd*JM,1+0*Jac); %% effectively, Bd*JM*Jac*JM'*Bd'

%% Set plotting interpolation operator
Nf = floor(1.2*N); [zf,wf]=zwuni(Nf); Jf=interp_mat(zf,z);

%% Set timestep size
[lam_max_est]=est_lam_cfl(U,V,Rx); ldt_max = 0.5;
lam_max_est=max(1e-5,lam_max_est);
dt=ldt_max/lam_max_est;
