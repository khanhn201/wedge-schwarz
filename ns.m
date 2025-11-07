clear all

hdr;    % 2-D SEM multi-element
close all;

N=10;


nu=1; alpha=1.e-0;

Re=1./nu;

%% Set ICs, problem parameters as function of N
[U,V,T,z,w,Dh,X,Y,Grr,Grs,Gss,Bl,Xr,Rx,Jac,Q,glo_num,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,BC_all,dA,interpdata_top]=set_sem_all(N);

[U_tip,V_tip,T_tip,z_tip,w_tip,Dh_tip,X_tip,Y_tip,Grr_tip,Grs_tip,Gss_tip,Bl_tip,Xr_tip,Rx_tip, ...
Jac_tip,Q_tip,glo_num_tip,Mu_tip,Mv_tip,Mp_tip,Mt_tip,ifnull_tip, ...
unxa_v_tip,unya_v_tip,BC_all_tip,dA_tip,interpdata_tip]=set_sem_all_tipv03(10);

##% Plot mesh
##E1 = size(X,2); E2 = size(X_tip,2); N1 = N+1;
##
## figure; hold on; axis equal;
## for e = 1:E1
##    %Extract patch for element e
##    Xe = squeeze(X(:,e,:));
##    Ye = squeeze(Y(:,e,:));
##    %Plot grid lines
##    plot(Xe, Ye, 'b-');           % lines along R
##    plot(Xe', Ye', 'b-');         % lines along S
##
##    % node_id = 0;
##    % for j = 1:N1; for i = 1:N1;
##    %   node_id = node_id +1;
##    %   text(Xe(i,j), Ye(i,j),num2str(node_id), 'fontsize',14);
##    % end;end
##    % xc = mean(Xe(:));
##    % yc = mean(Ye(:));
##
##    % add element number
##    % text(xc, yc, num2str(e), ...
##    %      'HorizontalAlignment','center', ...
##    %      'VerticalAlignment','middle', ...
##    %      'FontWeight','bold', ...
##    %      'Color','r');
##
## end
##
## for e = 1:E2
##    %Extract patch for element e
##    Xe = squeeze(X_tip(:,e,:));
##    Ye = squeeze(Y_tip(:,e,:));
##    %Plot grid lines
##    plot(Xe, Ye, 'r-');           % lines along R
##    plot(Xe', Ye', 'r-');         % lines along S
##
##    % node_id = 0;
##    % for j = 1:N1; for i = 1:N1;
##    %   node_id = node_id +1;
##    %   text(Xe(i,j), Ye(i,j),num2str(node_id), 'fontsize',14);
##    % end;end
##    % xc = mean(Xe(:));
##    % yc = mean(Ye(:));
##    %
##    % % add element number
##    % text(xc, yc, num2str(e), ...
##    %      'HorizontalAlignment','center', ...
##    %      'VerticalAlignment','middle', ...
##    %      'FontWeight','bold', ...
##    %      'Color','r');
##
## end
## pause;
%% Set dealiasing operators, JM,DM,BMh
[JM,DM,BMh,Jf,dt] = set_dealiasing(N,z,Jac,U,V,Rx);
[JM_tip,DM_tip,BMh_tip,Jf_tip,dt] = set_dealiasing(N,z_tip,Jac_tip,U_tip,V_tip,Rx_tip);




Tfinal = 4*pi; nsteps = ceil(Tfinal/dt)
dt = Tfinal/nsteps;
dt=1e-3; nsteps=999;

%% Initialize BDFk/EXTk arrays

P=0*X;  %% Initialize pressure to zero
P_tip= 0*X_tip;

%% Initialize BDFk/EXTk arrays
% O=0*X; F=O; G=O; H=O;
% U1=O;U2=O;U3=O; V1=O;V2=O;V3=O;
% F1=O;F2=O;F3=O; G1=O;G2=O;G3=O;
% f1=O;f2=O;f3=O; g1=O;g2=O;g3=O;
% T1=O;T2=O;T3=O; H1=O;H2=O;H3=O;

%% System-solve parameters
ifnull=0; tol=1.e-6; max_iter=140;

%%%%% TIME STEPPING LOOP %%%%%%

kk=0; k=0; time=0;
for iloop=1:1;
  figure;
  for istep =1:nsteps; k=k+1;

        time = time + dt;
        %% Interpolation nodes

        x_interp_top = interpdata_top(:,4);
        y_interp_top = interpdata_top(:,5);
        [Uinterp_top,Vinterp_top,Pinterp_top,Tinterp_top] = interpolate(x_interp_top,y_interp_top,X_tip,Y_tip,U_tip,V_tip,P_tip,T_tip,z_tip);

        %Solve 1 timestep for top
        [U,V,P,T] = solve_2dnse(N,U,V,P,T,Dh,X,Y,Grr,Grs,Gss,Bl,Rx,Jac,Q,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,dA,dt,JM,DM,BMh,istep,nu,alpha, ...
                               Uinterp_top,Vinterp_top,Pinterp_top,Tinterp_top,interpdata_top);

        x_interp_tip = interpdata_tip(:,4);
        y_interp_tip = interpdata_tip(:,5);
        [Uinterp_tip,Vinterp_tip,Pinterp_tip,Tinterp_tip] = interpolate(x_interp_tip,y_interp_tip,X,Y,U,V,P,T,z);

##        if mod(istep,1)==0 || istep==1;  kk=kk+1;
##
##    %      disp([itp itu itv itt])
##
##           hold off;
##
##           umax = max(max(max(abs(U))));
##           vmax = max(max(max(abs(V))));
##           tmax = max(max(max(abs(T))));
##           um(kk) = umax; vm(kk) = vmax;
##           tm(kk) = tmax; ti(kk) = time;
##
##           Uf = tensor3(Jf,1,Jf,U);  Uf=U;
##           Vf = tensor3(Jf,1,Jf,V);  Vf=V;
##           Xf = tensor3(Jf,1,Jf,X);  Xf=X;
##           Yf = tensor3(Jf,1,Jf,Y);  Yf=Y;
##           Tf = tensor3(Jf,1,Jf,T);  Tf=P;
##
##    %      if istep>2; Tf=Tf-Tfl; end;
##    %      Tfl = tensor3(Jf,1,Jf,T);
##    %      tmax = max(max(max(abs(Tf))));
##
##           s=['Time,UVT_{max}: ' num2str(time) ',   ' num2str(tmax) ,...
##              ', ' num2str(istep)'.'];
##    ##        se_mesh  (Xf,Yf,Tf,s);  hold on;
##           hold off; se_quiver(Xf,Yf,Uf,Vf,s);  axis equal; hold on;
##           % drawnow
##           time
##           se_mesh(X,Y,T, s)
##           drawnow; pause(0.1);
##
##    ##       hold on; se_quiver(Xf,Yf,U3plt,V3plt,s); drawnow;pause(0.1);
##    %      disp([umax vmax tmax])
##
##       end;


        %Solve for a tip
        [U_tip,V_tip,P_tip,T_tip] = solve_2dnse_tip(N,U_tip,V_tip,P_tip,T_tip,Dh_tip,X_tip,Y_tip,Grr_tip,Grs_tip,Gss_tip,Bl_tip,Rx_tip,Jac_tip,Q_tip,...
                                               Mu_tip,Mv_tip,Mp_tip,Mt_tip,ifnull_tip,unxa_v_tip,unya_v_tip,dA_tip,dt,JM_tip,DM_tip,BMh_tip, ...
                                               istep, nu, alpha,Uinterp_tip,Vinterp_tip,Pinterp_tip,Tinterp_tip,interpdata_tip);


  %    Diagonostics
##        U = U_tip; V = V_tip; P = P_tip; T = T_tip;
##        X = X_tip; Y = Y_tip; Jf = Jf_tip;

       if mod(istep,1)==0 || istep==1;  kk=kk+1;

  %      disp([itp itu itv itt])

         hold off;

  %      if istep>2; Tf=Tf-Tfl; end;
  %      Tfl = tensor3(Jf,1,Jf,T);
  %      tmax = max(max(max(abs(Tf))));

         % s=['Time,UVT_{max}: ' num2str(time) ',   ' num2str(tmax) ,...
         %    ', ' num2str(istep)'.'];
         s='Time'
         % se_mesh  (X,Y,P,s);hold on;
         % se_mesh  (X_tip,Y_tip,U_tip,s);
         % hold off; se_quiver(X,Y,U,V,s);  axis equal; hold on;
         % se_quiver(X_tip,Y_tip,U_tip,V_tip,s);  axis equal;
         x_all = [];
         y_all = [];
         u_all = [];
         v_all = [];
         N_top = size(X,1);
         for e = 1:size(X,2)
             x = reshape(X(:,e,:), N_top, N_top);
             y = reshape(Y(:,e,:), N_top, N_top);
             u = reshape(U(:,e,:), N_top, N_top);
             v = reshape(V(:,e,:), N_top, N_top);
             x_all = [x_all; x(:)];
             y_all = [y_all; y(:)];
             u_all = [u_all; u(:)];
             v_all = [v_all; v(:)];
         end
         N_tip = size(X_tip,1);
         for e = 1:size(X_tip,2)
             x = reshape(X_tip(:,e,:), N_tip, N_tip);
             y = reshape(Y_tip(:,e,:), N_tip, N_tip);
             u = reshape(U_tip(:,e,:), N_tip, N_tip);
             v = reshape(V_tip(:,e,:), N_tip, N_tip);
             x_all = [x_all; x(:)];
             y_all = [y_all; y(:)];
             u_all = [u_all; u(:)];
             v_all = [v_all; v(:)];
         end
         mag = sqrt(u_all.^2 + v_all.^2);
         mag(mag == 0) = 1;  % avoid division by zero
         u_all = u_all ./ mag;
         v_all = v_all ./ mag;
         % scatter(x_all, y_all, 10, mag, 'filled'); hold on;
         % colormap(jet);
         quiver(x_all, y_all, u_all, v_all, 'k');
         xlim([-0.5,0.5]);
         ylim([-0.1,1.2]);


         drawnow
         time


  ##       hold on; se_quiver(Xf,Yf,U3plt,V3plt,s); drawnow;pause(0.1);
  %      disp([umax vmax tmax])

       end;
##
##       if tmax > 10.9; break; end;
##       if umax > 10.9; break; end;
##       if vmax > 10.9; break; end;

   end;
end;

%%%%% TIME STEPPING LOOP %%%%%%

% kk=0; k=0; time=0;
% for iloop=1:1;
%   for istep =1:nsteps; k=k+1;
%
%      time = time+dt;
%
% %%   Set updated BDFk/EXTk coefficients
%      ndt = nu*dt; adt = alpha*dt;
%      if k==1; a1=1; a2=0; a3=0; b0=1; b1=1; b2=0; b3=0; end;
%      if k==2; a1=1.5; a2=-.5; a3=0; b0=1.5; b1=2; b2=-.5; b3=0; end;
%      if k==3; a1=3; a2=-3; a3=1; b0=11/6; b1=3; b2=-1.5; b3=2/6; end;
%      d1=dt*a1; d2=dt*a2; d3=dt*a3;
%
% %%   Set dealiased advecting field
%      [Cr,Cs]=set_advect_c(U,V,JM,BMh,Jac,Rx);
%
% %%   Set body force, volumetric heating
%      if k==1; QT =  0*Y; end;  %% No forcing, no heating
%      if k==1; FX =  0*Y; end;
%      if k==1; FY =  0*Y; end;
%
% %%   Evaluate curl-curl term (to be extrapolated)
%    [curlcurlX,curlcurlY,Omega]=curlcurl(U,V,Bl,Rx,Dh);
%      % curlcurlX = 0*X;
%      % curlcurlY = 0*Y;
% %    Omega = Lxi*(Dhx*V) - Lyi*(U*Dhy');
% %    curlcurlX =  Bl.*(Lyi*(Omega*Dhy'));
% %    curlcurlY = -Bl.*(Lxi*(Dhx*Omega));
%
% %
% %    Set Dirichlet conditions onto old fields
% %
%      [Ub,Vb,Tb]=set_dirichlet(U,V,T,Mu,Mv,Mt,X,Y);
%
%
% %%   Compute u-hat and u-tilde
%      U3=U2;U2=U1;U1=U;
%         F3=F2;F2=F1;F1=-advectl(U,Cr,Cs,JM,DM)+Bl.*FX;
%         f3=f2;f2=f1;f1=-nu*curlcurlX;
%         Uh=Bl.*(b1*U1+b2*U2+b3*U3)+(d1*F1+d2*F2+d3*F3);
%         Ut=Uh+(d1*f1+d2*f2+d3*f3);
%         Uh=Uh-axl(Ub,b0,ndt,Bl,Grr,Grs,Gss,Dh);
%
% %%   Compute v-hat and v-tilde
%      V3=V2;V2=V1;V1=V;
%         G3=G2;G2=G1;G1=-advectl(V,Cr,Cs,JM,DM)+Bl.*FY;
%         g3=g2;g2=g1;g1=-nu*curlcurlY;
%         Vh=Bl.*(b1*V1+b2*V2+b3*V3)+(d1*G1+d2*G2+d3*G3);
%         Vt=Vh+(d1*g1+d2*g2+d3*g3);
%         Vh=Vh-axl(Vb,b0,ndt,Bl,Grr,Grs,Gss,Dh);
%
% %%   Compute t-hat
%      T3=T2;T2=T1;T1=T;
%         H3=H2;H2=H1;H1=-advectl(T,Cr,Cs,JM,DM)+Bl.*QT;
%         Th=Bl.*(b1*T1+b2*T2+b3*T3)+(d1*H1+d2*H2+d3*H3);
%         Th=Th-axl(Tb,b0,adt,Bl,Grr,Grs,Gss,Dh);
%
% %    Pressure correction
%
%      divUt = weak_div(Ut,Vt,1.,Rx,Dh)/dt;
% %%   Add inhomogeneous Neumann data to divUT, if any. (Eq.(15) in split_slides.pdf)
%      b0dt = b0/dt;
%      divUt = divUt - b0dt*( (1-Mu).*unxa_v.*Ub - (1-Mv).*unya_v.*Vb );
%
% %%   Pressure-Poisson solve
%      h1=1; h0=0;
%      divUt = divUt -axl(P,h0,h1,Bl,Grr,Grs,Gss,Dh);
%      [dP,itp,res,lamda_h]=...
%          pcg_lambda(divUt,tol,max_iter,h0,h1,Mp,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull);
%      s=['Pressure. Step/Iter: = ' int2str([istep itp])];
% %    hold off; se_mesh  (X,Y,dP,s);  drawnow;
%      P = P+dP;
%
%      [dPdx,dPdy]=grad(P,Rx,Dh);
%      Uh = Uh - dt*Bl.*dPdx;
%      Vh = Vh - dt*Bl.*dPdy;
%
% %    Viscous/diffusive solves (diagonally-preconditioned CG):
%
% %%   Implicit solve - diagonal preconditioner
%      dAT=1./dA; dAT=1./(b0*qqt(Q,Bl)+adt*dAT);
%      dAU=1./dA; dAU=1./(b0*qqt(Q,Bl)+adt*dAU);
%      [U,itu,res,lamda_h]=...
%         pcg_lambda(Uh,tol,max_iter,b0,ndt,Mu,Q,Bl,Grr,Grs,Gss,Dh,dAU,ifnull);
%      [V,itv,res,lamda_h]=...
%         pcg_lambda(Vh,tol,max_iter,b0,ndt,Mv,Q,Bl,Grr,Grs,Gss,Dh,dAU,ifnull);
%      [T,itt,res,lamda_h]=...
%         pcg_lambda(Th,tol,max_iter,b0,adt,Mt,Q,Bl,Grr,Grs,Gss,Dh,dAT,ifnull);
%
%      U=U+Ub;  %% Add back any prescribed Dirichlet conditions
%      V=V+Vb;
%      T=T+Tb;
%
% %    Diagonostics
%      if mod(istep,5)==0 || istep==1;  kk=kk+1;
%        hold off;
%        se_mesh  (X,Y,U,'s');hold on;
%
% %      disp([itp itu itv itt])
%
% %
% %        umax = max(max(max(abs(U))));
% %        vmax = max(max(max(abs(V))));
% %        tmax = max(max(max(abs(T))));
% %        um(kk) = umax; vm(kk) = vmax;
% %        tm(kk) = tmax; ti(kk) = time;
% %
% %        Uf = tensor3(Jf,1,Jf,U);  Uf=U;
% %        Vf = tensor3(Jf,1,Jf,V);  Vf=V;
% %        Xf = tensor3(Jf,1,Jf,X);  Xf=X;
% %        Yf = tensor3(Jf,1,Jf,Y);  Yf=Y;
% %        Tf = tensor3(Jf,1,Jf,T);  Tf=P;
% %
% % %      if istep>2; Tf=Tf-Tfl; end;
% % %      Tfl = tensor3(Jf,1,Jf,T);
% % %      tmax = max(max(max(abs(Tf))));
% %
% %        s=['Time,UVT_{max}: ' num2str(time) ',   ' num2str(tmax) ,...
% %           ', ' num2str(istep) ', ' num2str(itp)'.'];
% %        % se_mesh  (Xf,Yf,Tf,s);  hold on;
% %        hold off; se_quiver(Xf,Yf,Uf,Vf,s);  axis equal; hold on;
% %        % drawnow
% %        time
% %        se_mesh(X,Y,T, s)
% %        drawnow
% %      disp([umax vmax tmax])
%
%      end;
%
%      if tmax > 1.9; break; end;
%      if umax > 1.9; break; end;
%      if vmax > 1.9; break; end;
%
%    end;
% end;
%
%
