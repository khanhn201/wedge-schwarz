function [Ufinal,Vfinal,Pfinal,Tfinal]= solve_2dnse_tip(N,U,V,P,T,Dh,X,Y,Grr,Grs,Gss,Bl,Rx,Jac,Q,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,dA,dt,JM,DM,BMh,istep, nu, alpha)

k = istep;


persistent U1 = 0*X; persistent U2 = 0*X; persistent U3 = 0*X;
persistent V1 = 0*X; persistent V2 = 0*X; persistent V3 = 0*X;

persistent F1 = 0*X; persistent F2 = 0*X; persistent F3 = 0*X;
persistent f1 = 0*X; persistent f2 = 0*X; persistent f3 = 0*X;


persistent G1 = 0*X; persistent G2 = 0*X; persistent G3 = 0*X;
persistent g1 = 0*X; persistent g2 = 0*X; persistent g3 = 0*X;


persistent T1 = 0*X; persistent T2 = 0*X; persistent T3 = 0*X;
persistent H1 = 0*X; persistent H2 = 0*X; persistent H3 = 0*X;


##F1=O;F2=O;F3=O; G1=O;G2=O;G3=O;
##f1=O;f2=O;f3=O; g1=O;g2=O;g3=O;
##T1=O;T2=O;T3=O; H1=O;H2=O;H3=O;
%% System-solve parameters
ifnull=0; tol=1.e-6; max_iter=140;

%%   Set updated BDFk/EXTk coefficients
     ndt = nu*dt; adt = alpha*dt;
     if k==1; a1=1; a2=0; a3=0; b0=1; b1=1; b2=0; b3=0; end;
     if k==2; a1=1.5; a2=-.5; a3=0; b0=1.5; b1=2; b2=-.5; b3=0; end;
     if k>=3; a1=3; a2=-3; a3=1; b0=11/6; b1=3; b2=-1.5; b3=2/6; end;
     d1=dt*a1; d2=dt*a2; d3=dt*a3;

%%   Set dealiased advecting field
     [Cr,Cs]=set_advect_c(U,V,JM,BMh,Jac,Rx);

%%   Set body force, volumetric heating
     if k>=1; QT =  0*Y; end;  %% No forcing, no heating
     if k>=1; FX =  0*Y; end;
     if k>=1; FY =  0*Y; end;

%%   Evaluate curl-curl term (to be extrapolated)
    [curlcurlX,curlcurlY,Omega]=curlcurl(U,V,Bl,Rx,Dh);
##    Omega = Lxi*(Dhx*V) - Lyi*(U*Dhy');
##    curlcurlX =  Bl.*(Lyi*(Omega*Dhy'));
##    curlcurlY = -Bl.*(Lxi*(Dhx*Omega));

%
%    Set Dirichlet conditions onto old fields
%
     [Ub,Vb,Tb]=set_dirichlet(U,V,T,Mu,Mv,Mt,X,Y);


%%   Compute u-hat and u-tilde
     U3=U2;U2=U1;U1=U;
        F3=F2;F2=F1;F1=-advectl(U,Cr,Cs,JM,DM)+Bl.*FX;
        f3=f2;f2=f1;f1=-nu*curlcurlX;
        Uh=Bl.*(b1*U1+b2*U2+b3*U3)+(d1*F1+d2*F2+d3*F3);
        Ut=Uh+(d1*f1+d2*f2+d3*f3);
        Uh=Uh-axl(Ub,b0,ndt,Bl,Grr,Grs,Gss,Dh);

%%   Compute v-hat and v-tilde
     V3=V2;V2=V1;V1=V;
        G3=G2;G2=G1;G1=-advectl(V,Cr,Cs,JM,DM)+Bl.*FY;
        g3=g2;g2=g1;g1=-nu*curlcurlY;
        Vh=Bl.*(b1*V1+b2*V2+b3*V3)+(d1*G1+d2*G2+d3*G3);
        Vt=Vh+(d1*g1+d2*g2+d3*g3);
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
     divUt = divUt - b0dt*( (1-Mu).*unxa_v.*Ub - (1-Mv).*unya_v.*Vb );

%%   Pressure-Poisson solve
     h1=1; h0=0;
     divUt = divUt -axl(P,h0,h1,Bl,Grr,Grs,Gss,Dh);
     [dP,itp,res,lamda_h]=...
         pcg_lambda(divUt,tol,max_iter,h0,h1,Mp,Q,Bl,Grr,Grs,Gss,Dh,dA,ifnull);
     s=['Pressure. Step/Iter: = ' int2str([istep itp])];
%    hold off; se_mesh  (X,Y,dP,s);  drawnow;
     Pfinal = P+dP;

     [dPdx,dPdy]=grad(P,Rx,Dh);
     Uh = Uh - dt*Bl.*dPdx;
     Vh = Vh - dt*Bl.*dPdy;

%    Viscous/diffusive solves (diagonally-preconditioned CG):

%%   Implicit solve - diagonal preconditioner
     dAT=1./dA; dAT=1./(b0*qqt(Q,Bl)+adt*dAT);
     dAU=1./dA; dAU=1./(b0*qqt(Q,Bl)+adt*dAU);
     [U,itu,res,lamda_h]=...
        pcg_lambda(Uh,tol,max_iter,b0,ndt,Mu,Q,Bl,Grr,Grs,Gss,Dh,dAU,ifnull);
     [V,itv,res,lamda_h]=...
        pcg_lambda(Vh,tol,max_iter,b0,ndt,Mv,Q,Bl,Grr,Grs,Gss,Dh,dAU,ifnull);
     [T,itt,res,lamda_h]=...
        pcg_lambda(Th,tol,max_iter,b0,adt,Mt,Q,Bl,Grr,Grs,Gss,Dh,dAT,ifnull);

     Ufinal=U+Ub;  %% Add back any prescribed Dirichlet conditions
     Vfinal=V+Vb;
     Tfinal=T+Tb;



