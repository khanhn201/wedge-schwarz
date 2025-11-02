
hdr;    % 2-D SEM multi-element
close all;

N=4;

nu=1.e-1/4; alpha=1.e-4;

Re=1./nu;

%% Set ICs, problem parameters as function of N
[U,V,T,z,w,Dh,X,Y,Grr,Grs,Gss,Bl,Xr,Rx,Jac,Q,glo_num,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,BC_all,dA]=set_sem_all(N);
[U_tip,V_tip,T_tip,z_tip,w_tip,Dh_tip,X_tip,Y_tip,Grr_tip,Grs_tip,Gss_tip,Bl_tip,Xr_tip,Rx_tip, ...
Jac_tip,Q_tip,glo_num_tip,Mu_tip,Mv_tip,Mp_tip,Mt_tip,ifnull_tip, ...
unxa_v_tip,unya_v_tip,BC_all_tip,dA_tip]=set_sem_all_tipv02(N);

##%% Set dealiasing operators, JM,DM,BMh
##Nd = floor(1.5*N);
##[zd,wd]=zwgl(Nd); Bd=diag(wd);
##JM=interp_mat(zd,z);
##DM=deriv_mat (zd);
##BMh=tensor3(Bd*JM,1,Bd*JM,1+0*Jac); %% effectively, Bd*JM*Jac*JM'*Bd'
##
##%% Set plotting interpolation operator
##Nf = floor(1.2*N); [zf,wf]=zwuni(Nf); Jf=interp_mat(zf,z);
##
##%% Set timestep size
##[lam_max_est]=est_lam_cfl(U,V,Rx); ldt_max = 0.5;
##lam_max_est=max(1e-5,lam_max_est);
##dt=ldt_max/lam_max_est;

[JM,DM,BMh,Jf,dt] = set_dealiasing(N,z,Jac,U,V,Rx);
[JM_tip,DM_tip,BMh_tip,Jf_tip,dt] = set_dealiasing(N,z_tip,Jac_tip,U_tip,V_tip,Rx_tip);



dt=.1;
Tfinal = 4*pi; nsteps = ceil(Tfinal/dt)
dt = Tfinal/nsteps;
dt=.01; nsteps=999;
##
##%% Initialize BDFk/EXTk arrays
O=0*X_tip; F=O; G=O; H=O;

P=O;  %% Initialize pressure to zero
P_tip= O;

%% System-solve parameters
ifnull=0; tol=1.e-6; max_iter=140;

%%%%% TIME STEPPING LOOP %%%%%%

kk=0; k=0; time=0;
for iloop=1:1;
  for istep =1:nsteps; k=k+1;

      time = dt*istep;
      %Solve 1 timestep for top
       [U,V,P,T,U3plt,V3plt] = solve_2dnse(N,U,V,P,T,Dh,X,Y,Grr,Grs,Gss,Bl,Rx,Jac,Q,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,dA,dt,JM,DM,BMh,istep,nu,alpha);

      %Solve for a tip
##      [U_tip,V_tip,P_tip,T_tip] = solve_2dnse_tip(N,U_tip,V_tip,P_tip,T_tip,Dh_tip,X_tip,Y_tip,Grr_tip,Grs_tip,Gss_tip,Bl_tip,Rx_tip,Jac_tip,Q_tip,...
##                                             Mu_tip,Mv_tip,Mp_tip,Mt_tip,ifnull_tip,unxa_v_tip,unya_v_tip,dA_tip,dt,JM_tip,DM_tip,BMh_tip, ...
##                                             istep, nu, alpha);

%    Diagonostics
##      U = U_tip; V = V_tip; P = P_tip; T = T_tip;
##      X = X_tip; Y = Y_tip; Jf = Jf_tip;
     if mod(istep,5)==0 || istep==1;  kk=kk+1;

%      disp([itp itu itv itt])

       hold off;

       umax = max(max(max(abs(U))));
       vmax = max(max(max(abs(V))));
       tmax = max(max(max(abs(T))));
       um(kk) = umax; vm(kk) = vmax;
       tm(kk) = tmax; ti(kk) = time;

       Uf = tensor3(Jf,1,Jf,U);  Uf=U;
       Vf = tensor3(Jf,1,Jf,V);  Vf=V;
       Xf = tensor3(Jf,1,Jf,X);  Xf=X;
       Yf = tensor3(Jf,1,Jf,Y);  Yf=Y;
       Tf = tensor3(Jf,1,Jf,T);  Tf=P;

%      if istep>2; Tf=Tf-Tfl; end;
%      Tfl = tensor3(Jf,1,Jf,T);
%      tmax = max(max(max(abs(Tf))));

       s=['Time,UVT_{max}: ' num2str(time) ',   ' num2str(tmax) ,...
          ', ' num2str(istep)'.'];
##        se_mesh  (Xf,Yf,Tf,s);  hold on;
       hold off; se_quiver(Xf,Yf,Uf,Vf,s);  axis equal; hold on;0
       % drawnow
       time
       se_mesh(X,Y,T, s)
       drawnow; pause(0.1);

##       hold on; se_quiver(Xf,Yf,U3plt,V3plt,s); drawnow;pause(0.1);
%      disp([umax vmax tmax])

     end;

     if tmax > 10.9; break; end;
     if umax > 10.9; break; end;
     if vmax > 10.9; break; end;

   end;
end;

