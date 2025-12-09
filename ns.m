clear all

hdr;    % 2-D SEM multi-element
close all;

N=6;


nu=1.0/5.0; alpha=1.e-0;

Re=1./nu;

%% Set ICs, problem parameters as function of N
[U,V,T,z,w,Dh,X,Y,Grr,Grs,Gss,Bl,Xr,Rx,Jac,Q,glo_num,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,BC_all,dA,interpdata_top]=set_sem_all(N);

[U_tip,V_tip,T_tip,z_tip,w_tip,Dh_tip,X_tip,Y_tip,Grr_tip,Grs_tip,Gss_tip,Bl_tip,Xr_tip,Rx_tip, ...
Jac_tip,Q_tip,glo_num_tip,Mu_tip,Mv_tip,Mp_tip,Mt_tip,ifnull_tip, ...
unxa_v_tip,unya_v_tip,BC_all_tip,dA_tip,interpdata_tip]=set_sem_all_tip(4);

omegx = -Y_tip; omegy = X_tip;
#Plot mesh
% E1 = size(X,2); E2 = size(X_tip,2); N1 = N+1;
% figure; hold on; axis equal;
% for e = 1:E1
%     %Extract patch for element e
%     Xe = squeeze(X(:,e,:));
%     Ye = squeeze(Y(:,e,:));
%     %Plot grid lines
%     plot(Xe, Ye, 'b-');           % lines along R
%     plot(Xe', Ye', 'b-');         % lines along S
%
%      node_id = 0;
%      % for j = 1:N1; for i = 1:N1;
%      %   node_id = node_id +1;
%      %   text(Xe(i,j), Ye(i,j),num2str(node_id), 'fontsize',14);
%      % end;end
%      xc = mean(Xe(:));
%      yc = mean(Ye(:));
%
% ##     add element number
%      % text(xc, yc, num2str(e), ...
%      %      'HorizontalAlignment','center', ...
%      %      'VerticalAlignment','middle', ...
%      %      'FontWeight','bold', ...
%      %      'Color','r');
%
%  end
%
%  for e = 1:E2
%     %Extract patch for element e
%     Xe = squeeze(X_tip(:,e,:));
%     Ye = squeeze(Y_tip(:,e,:));
%     %Plot grid lines
%     plot(Xe, Ye, 'r-');           % lines along R
%     plot(Xe', Ye', 'r-');         % lines along S
%
%      node_id = 0;
%      % for j = 1:N1; for i = 1:N1;
%      %   node_id = node_id +1;
%      %   text(Xe(i,j), Ye(i,j),num2str(node_id), 'fontsize',14);
%      % end;end
%      xc = mean(Xe(:));
%      yc = mean(Ye(:));
%
%      % add element number
%      % text(xc, yc, num2str(e), ...
%      %      'HorizontalAlignment','center', ...
%      %      'VerticalAlignment','middle', ...
%      %      'FontWeight','bold', ...
%      %      'Color','r');
%
%  end
%  pause;
%% Set dealiasing operators, JM,DM,BMh
[JM,DM,BMh,Jf,dt] = set_dealiasing(N,z,Jac,U,V,Rx);
[JM_tip,DM_tip,BMh_tip,Jf_tip,dt] = set_dealiasing(N,z_tip,Jac_tip,U_tip,V_tip,Rx_tip);




Tfinal = 4*pi; nsteps = ceil(Tfinal/dt)
dt = Tfinal/nsteps;
dt=1e-1; nsteps=999;

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

        time = time + dt
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


        %Solve for a tip
        for k = 1:1
        [U_tip,V_tip,P_tip,T_tip] = solve_2dnse_tip(N,U_tip,V_tip,P_tip,T_tip,Dh_tip,X_tip,Y_tip,Grr_tip,Grs_tip,Gss_tip,Bl_tip,Rx_tip,Jac_tip,Q_tip,...
                                               Mu_tip,Mv_tip,Mp_tip,Mt_tip,ifnull_tip,unxa_v_tip,unya_v_tip,dA_tip,dt/1.0,JM_tip,DM_tip,BMh_tip, ...
                                               istep, nu, alpha,Uinterp_tip,Vinterp_tip,Pinterp_tip,Tinterp_tip,interpdata_tip, omegx, omegy);
        end

        c = cos(dt);
        s = sin(dt);
  
        Xnew =  c*X_tip - s*Y_tip;
        Ynew =  s*X_tip + c*Y_tip;
        X_tip = Xnew;
        Y_tip = Ynew;

        Xintnew = c*interpdata_tip(:,4) - s*interpdata_tip(:,5);
        Yintnew = s*interpdata_tip(:,4) + c*interpdata_tip(:,5);
        interpdata_tip(:,4) = Xintnew;
        interpdata_tip(:,5) = Yintnew;

        omegx = -Y_tip;
        omegy = X_tip;
  %    Diagonostics
##        U = U_tip; V = V_tip; P = P_tip; T = T_tip;
##        X = X_tip; Y = Y_tip; Jf = Jf_tip;

       #if mod(istep,1)==0 || istep==1;  kk=kk+1;
        if (1 ==1);

  %      disp([itp itu itv itt])
        [Fx_total, Fy_total] = compute_drag_lift(U_tip,V_tip,P_tip,nu,Xr_tip,Rx_tip,Dh_tip,w_tip);
        Fy_total

         hold off;

  %      if istep>2; Tf=Tf-Tfl; end;
  %      Tfl = tensor3(Jf,1,Jf,T);
  %      tmax = max(max(max(abs(Tf))));

         % s=['Time,UVT_{max}: ' num2str(time) ',   ' num2str(tmax) ,...
         %    ', ' num2str(istep)'.'];
##          se_mesh  (X,Y,P,s);hold on;
##          se_mesh  (X_tip,Y_tip,P_tip,s);zlim([-0.3,0.3])
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
##         u_all = u_all ./ mag;
##         v_all = v_all ./ mag;
         scatter(x_all, y_all, 10, mag, 'filled'); hold on;
         colormap(jet);
        quiver(x_all, y_all, u_all, v_all, 'k');
        drawnow;




  ##       hold on; se_quiver(Xf,Yf,U3plt,V3plt,s); drawnow;pause(0.1);
  %      disp([umax vmax tmax])

       end;
##
##       if tmax > 10.9; break; end;
##       if umax > 10.9; break; end;
##       if vmax > 10.9; break; end;

   end;
end;


