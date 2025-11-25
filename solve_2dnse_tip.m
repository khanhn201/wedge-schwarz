function [U,V,P,T]= solve_2dnse_tip(N,U,V,P,T,Dh,X,Y,Grr,Grs,Gss,Bl,Rx,Jac,Q,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,dA,dt,JM,DM,BMh,istep, nu, alpha, ...
                                                        Uinterp_tip,Vinterp_tip,Pinterp_tip,Tinterp_tip,interpdata_tip)


%%[U,V,P,T,U3plt,V3plt] = solve_2dnse(N,U,V,P,T,Dh,X,Y,Grr,Grs,Gss,Bl,Rx,Jac,Q,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,dA,dt,JM,DM,BMh,istep,nu,alpha);
%k = istep;
k = 1;
N1 = size(X,1);
E = size(X,2);
%% System-solve parameters
tol=1.e-16; max_iter=5000;


%%persistent U1 = 0*X; persistent U2 = 0*X; persistent U3 = 0*X;
%%persistent V1 = 0*X; persistent V2 = 0*X; persistent V3 = 0*X;
%%persistent F1 = 0*X; persistent F2 = 0*X; persistent F3 = 0*X;
%%persistent f1 = 0*X; persistent f2 = 0*X; persistent f3 = 0*X;
%%persistent G1 = 0*X; persistent G2 = 0*X; persistent G3 = 0*X;
%%persistent g1 = 0*X; persistent g2 = 0*X; persistent g3 = 0*X;
%%persistent T1 = 0*X; persistent T2 = 0*X; persistent T3 = 0*X;
%%persistent H1 = 0*X; persistent H2 = 0*X; persistent H3 = 0*X;
%% --- Persistent history buffers (INIT SAFELY) ---
persistent U1 U2 U3 V1 V2 V3 T1 T2 T3
persistent F1 F2 F3 G1 G2 G3 H1 H2 H3
persistent f1 f2 f3 g1 g2 g3           % (viscous/curlcurl parts)
persistent last_sz                        % to detect grid changes
persistent P_prev

% Bootstrap or re-bootstrap if:
%  - first call (isempty)
%  - istep == 1 (new time integration)
%  - grid size changed
need_init = isempty(U1) || k==1 || isempty(last_sz) || ~isequal(last_sz, size(U));

if need_init
    Z  = 0*U;     % zero of correct size
    U1 = Z; U2 = Z; U3 = Z;
    V1 = Z; V2 = Z; V3 = Z;
    T1 = Z; T2 = Z; T3 = Z;

    F1 = Z; F2 = Z; F3 = Z;    % conv/forcing history for U
    G1 = Z; G2 = Z; G3 = Z;    % conv/forcing history for V
    H1 = Z; H2 = Z; H3 = Z;    % conv/forcing history for T

    f1 = Z; f2 = Z; f3 = Z;    % viscous (curlcurlX) history for U
    g1 = Z; g2 = Z; g3 = Z;    % viscous (curlcurlY) history for V

    last_sz = size(U);
end


  %%   Set updated BDFk/EXTk coefficients
    nu
    alpha
    dt
     ndt = nu*dt; adt = alpha*dt;
     if k==1; a1=1; a2=0; a3=0; b0=1; b1=1; b2=0; b3=0; end;
     if k==2; a1=1.5; a2=-.5; a3=0; b0=1.5; b1=2; b2=-.5; b3=0; end;
     if k>=3; a1=3; a2=-3; a3=1; b0=11/6; b1=3; b2=-1.5; b3=2/6; end;
     d1=dt*a1; d2=dt*a2; d3=dt*a3;


%%     printf("-------- solve_2dnse ---------- \n");
%%     printf("k = %f \n",k);
%%     printf("a1= %f, a2 = %f, a3 = %f | b0 = %f, b1 = %f, b2 = %f, b3 = %f",a1, a2,a3,b0, b1,b2,b3);

%%   Set dealiased advecting field
     [Cr,Cs]=set_advect_c(U,V,JM,BMh,Jac,Rx);

%%   Set body force, volumetric heating
     if k>=1; QT =  0*Y; end;  %% No forcing, no heating
     if k>=1; FX =  0*Y; end;
     if k>=1; FY =  0*Y; end;


%%   Evaluate curl-curl term (to be extrapolated)
    [curlcurlX,curlcurlY,Omega]=curlcurl(U,V,Bl,Rx,Dh);
     # curlcurlX = 0*X;
     # curlcurlY = 0*Y;
%%    Omega = Lxi*(Dhx*V) - Lyi*(U*Dhy');
%%    curlcurlX =  Bl.*(Lyi*(Omega*Dhy'));
%%    curlcurlY = -Bl.*(Lxi*(Dhx*Omega));

%
%    Set Dirichlet conditions onto old fields
%
     [Ub,Vb,Pb,Tb]=set_dirichlet_tip(U,V,P,T,Mu,Mv,Mt,X,Y); %for tip part

     % interaction??
     e = interpdata_tip(:,1); r = interpdata_tip(:,2); s = interpdata_tip(:,3);
     for i = 1:size(interpdata_tip,1);
       Ub(r(i),e(i),s(i)) = Uinterp_tip(i);
       Vb(r(i),e(i),s(i)) = Vinterp_tip(i);
       Pb(r(i),e(i),s(i)) = Pinterp_tip(i);
       Tb(r(i),e(i),s(i)) = Tinterp_tip(i);
     endfor

%%   Compute u-hat and u-tilde
     U3=U2;U2=U1;U1=U;
        F3=F2;F2=F1;F1=-advectl(U,Cr,Cs,JM,DM)+Bl.*FX;
        f3=f2;f2=f1;f1=-nu*curlcurlX;
        Uh=Bl.*(b1*U1+b2*U2+b3*U3)+(d1*F1+d2*F2+d3*F3);
        Ut=Uh+(d1*f1+d2*f2+d3*f3);
        % Ut=Uh + f1;
        Uh=Uh-axl(Ub,b0,ndt,Bl,Grr,Grs,Gss,Dh);

%%   Compute v-hat and v-tilde
     V3=V2;V2=V1;V1=V;
        G3=G2;G2=G1;G1=-advectl(V,Cr,Cs,JM,DM)+Bl.*FY;
        g3=g2;g2=g1;g1=-nu*curlcurlY;
        Vh=Bl.*(b1*V1+b2*V2+b3*V3)+(d1*G1+d2*G2+d3*G3);
        Vt=Vh+(d1*g1+d2*g2+d3*g3);
        % Vt=Vh+g1;
        Vh=Vh-axl(Vb,b0,ndt,Bl,Grr,Grs,Gss,Dh);

%%   Compute t-hat
     T3=T2;T2=T1;T1=T;
        H3=H2;H2=H1;H1=-advectl(T,Cr,Cs,JM,DM)+Bl.*QT;
        Th=Bl.*(b1*T1+b2*T2+b3*T3)+(d1*H1+d2*H2+d3*H3);
        Th=Th-axl(Tb,b0,adt,Bl,Grr,Grs,Gss,Dh);

%    Pressure correction

     divUt = weak_div(Ut,Vt,1.,Rx,Dh)/dt;
%%   Add inhomogeneous Neumann data to divUT, if any. (Eq.(15) in split_slides.pdf)
     b0dt = b0/dt;
     divUt = divUt - b0dt*( (1-Mu).*unxa_v.*Ub + (1-Mv).*unya_v.*Vb );

%%   Pressure-Poisson solve
     h1=1; h0=0;
      divUt = divUt -axl(Pb,h0,h1,Bl,Grr,Grs,Gss,Dh);
      [dP,itp,res,lamda_h]=...
          pcg_lambda(divUt,tol,max_iter,h0,h1,Mp,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull);
      res
      s=['Pressure. Step/Iter: = ' int2str([istep itp])];
 %    hold off; se_mesh  (X,Y,dP,s);  drawnow;
      P = Pb+dP;



##     P_bar = 0*P;
##     for i=1:istep-1
##         Pp = reshape(P_prev(i,:, : , :), [N1 E N1]);
##         alphai = sum(sum(sum(Pp.*divUt)));
##         P_bar = P_bar + alphai*Pp;
##     end
##     divUt = divUt - axl(P_bar,h0,h1,Bl,Grr,Grs,Gss,Dh);
##     [dP,itp,res,lamda_h]=...
##         pcg_lambda(divUt,tol,max_iter,h0,h1,Mp,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull);
##     % s=['Pressure. Step/Iter: = ' int2str([istep itp])];
##     res
##     P = P_bar+dP;


     [dPdx,dPdy]=grad(P,Rx,Dh);
     Uh = Uh - dt*Bl.*dPdx;
     Vh = Vh - dt*Bl.*dPdy;

%    Viscous/diffusive solves (diagonally-preconditioned CG):

%%   Implicit solve - diagonal preconditioner
     dAT=1./dA; dAT=1./(b0*qqt(Q,Bl)+adt*dAT);
     dAU=1./dA; dAU=1./(b0*qqt(Q,Bl)+adt*dAU);
     [U,itu,res,lamda_h]=...
        pcg_lambda(Uh,tol,max_iter,b0,ndt,Mu,Q,Bl,Grr,Grs,Gss,Dh,dAU,ifnull);
     res
     [V,itv,res,lamda_h]=...
        pcg_lambda(Vh,tol,max_iter,b0,ndt,Mv,Q,Bl,Grr,Grs,Gss,Dh,dAU,ifnull);
     res
     [T,itt,res,lamda_h]=...
        pcg_lambda(Th,tol,max_iter,b0,adt,Mt,Q,Bl,Grr,Grs,Gss,Dh,dAT,ifnull);

     U=U+Ub;  %% Add back any prescribed Dirichlet conditions
     V=V+Vb;
     T=T+Tb;
##
##     P_tilde = P;
##     for i=1:istep-1
##         Pp = reshape(P_prev(i,:, : , :), [N1 E N1]);
##         alphai = sum(sum(sum(Pp.*axl(P,h0,h1,Bl,Grr,Grs,Gss,Dh))));
##
##         P_tilde = P_tilde - alphai*Pp;
##     end
##     beta = sum(sum(sum(P_tilde.*axl(P_tilde,h0,h1,Bl,Grr,Grs,Gss,Dh))));
##     if beta != 0
##         P_tilde =  P_tilde/beta;
##     end
##     P_prev(istep, :, :, :) = P_tilde;
