function [U,V,T,z,w,Dh,X,Y,Grr,Grs,Gss,Bl,Xr,Rx,Jac,Q,glo_num,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,BC_all,dA]...
             = set_sem_all_tipv02(N);

hdr;    % 2-D SEM multi-element

Nelx = 2;  Nely = 2; E = Nelx*Nely;
% Nelx = 1;  Nely = 1; E = Nelx*Nely;
N1=N+1;

%% Problem parameters, as function of N

N1=N+1;
##E=Nelx*Nely;

[z,w]=zwgll(N);                % Set basic operators
Dh=deriv_mat(z);

E = 1;
[R,S]=ndgrid(z,z);             % Build SEM mesh
X=zeros(N1,E,N1); Y=X;
e=0;
% geometric control
alpha = 90/2 * pi/180;%28.5/2 * pi/180;         % half-aperture angle in radians
y_top = 1.0;                   % top of trapezoid
y_mid = 0.0;                   % interface between trapezoid and kite
y_bot = -1.0;                  % apex (bottom)
bend  = 0.2;                   % inward bending of top centerline (optional)

% slopes from apex (constant on both sides)
half_width = @(y) abs(y - y_bot) * tan(alpha);

% -------------------------------------------------------------------------
##% === (1) Bent trapezoid: two halves merged together ===
##% Compute key widths from slope
##w_top = 2 * half_width(y_top);   % total width at top
##w_mid = 2 * half_width(y_mid);   % total width at interface
##
##% Left trapezoid
##x1L = -w_mid/2;  y1L = y_mid;
##x4L = -w_top/2;  y4L = y_top - 0.25;    % slightly shorter top
##x3L =  0*bend;   y3L = y_top;           % inward bent top
##x2L =  0.0;      y2L = y_mid + 0.5;     % joint centerline
##
##X_L = 0.25*((1-R).*(1-S)*x1L + (1+R).*(1-S)*x2L + ...
##            (1+R).*(1+S)*x3L + (1-R).*(1+S)*x4L);
##Y_L = 0.25*((1-R).*(1-S)*y1L + (1+R).*(1-S)*y2L + ...
##            (1+R).*(1+S)*y3L + (1-R).*(1+S)*y4L);
##
##% Right trapezoid (mirror)
##x1R =  w_mid/2;  y1R = y_mid;
##x2R =  w_top/2;  y2R = y_top - 0.25;
##x3R = 0*-bend;     y3R = y_top;
##x4R =  0.0;      y4R = y_mid + 0.5;
##
##X_R = 0.25*((1-R).*(1-S)*x1R + (1+R).*(1-S)*x2R + ...
##            (1+R).*(1+S)*x3R + (1-R).*(1+S)*x4R);
##Y_R = 0.25*((1-R).*(1-S)*y1R + (1+R).*(1-S)*y2R + ...
##            (1+R).*(1+S)*y3R + (1-R).*(1+S)*y4R);

% -------------------------------------------------------------------------
% === (2) Kite / apex element ===
% Ensure edges follow same slope
x1K = 0.0;        y1K = y_bot;       % bottom tip (apex)
x2K = half_width(y_mid);   y2K = y_mid;       % right joint
x3K = 0.0;        y3K = y_mid + 0.5; % upper center of kite
x4K = -half_width(y_mid);  y4K = y_mid;       % left joint

X_K = 0.25*((1-R).*(1-S)*x1K + (1+R).*(1-S)*x2K + ...
            (1+R).*(1+S)*x3K + (1-R).*(1+S)*x4K);
Y_K = 0.25*((1-R).*(1-S)*y1K + (1+R).*(1-S)*y2K + ...
            (1+R).*(1+S)*y3K + (1-R).*(1+S)*y4K);
##Y_K(:) = Y_K(:)*2+1;
% -------------------------------------------------------------------------
% Plot everything
% -------------------------------------------------------------------------
##figure; hold on; axis equal; box on;
##
##% top trapezoid (left + right)
##plot(X_L, Y_L, 'k'); plot(X_L', Y_L', 'k');
##plot(X_R, Y_R, 'k'); plot(X_R', Y_R', 'k');
##
##% bottom kite
##plot(X_K, Y_K, 'k'); plot(X_K', Y_K', 'k');
##
##% x_test= l
##
##xlabel('x'); ylabel('y');
##title(sprintf('Wedge Mesh with Aperture α = %.1f°', 28.5));

X(:,1,:)= X_K; Y(:,1,:) = Y_K;



 figure; hold on; axis equal;
 for e = 1:1
     % Extract patch for element e
     Xe = squeeze(X(:,e,:));
     Ye = squeeze(Y(:,e,:));
     % Plot grid lines
     plot(Xe, Ye, 'k-');           % lines along R
     plot(Xe', Ye', 'k-');         % lines along S

     node_id = 0;
     for j = 1:N1; for i = 1:N1;
       node_id = node_id +1;
       text(Xe(i,j), Ye(i,j),num2str(node_id), 'fontsize',14);
     end;end
     xc = mean(Xe(:));
     yc = mean(Ye(:));

     % add element number
     text(xc, yc, num2str(e), ...
          'HorizontalAlignment','center', ...
          'VerticalAlignment','middle', ...
          'FontWeight','bold', ...
          'Color','r');

 end
 pause;

 % [X,Y]=morph_circ(X,Y);         % Morph mesh

[Grr,Grs,Gss,Bl,Xr,Rx,Jac]=geom_elem(X,Y,Dh,w); % Terms for "A"
vol = sum(sum(sum(Bl)))
[Q,glo_num]=set_tp_semq_tip(N);


BC_all = [ 'D' 'N' 'D' 'D' ;     %% U
           'D' 'N' 'D' 'D' ;     %% V
           'N' 'D' 'N' 'N' ;     %% P
           'D' 'N' 'N' 'N' ];    %% T
Nelx = 1;Nely = 1;
[Mu,Q,glo_num]=set_mask(BC_all(1,:),1,1,Q,glo_num);
[Mv,Q,glo_num]=set_mask(BC_all(2,:),Nelx,Nely,Q,glo_num);
[Mp,Q,glo_num]=set_mask(BC_all(3,:),Nelx,Nely,Q,glo_num); ifnull=1;
[Mt,Q,glo_num]=set_mask(BC_all(4,:),Nelx,Nely,Q,glo_num);

[unxa_v,unya_v] = set_unxy(Mu,X,Y,Xr);

dA=diag_sem(Grr,Grs,Gss,Dh); dA=qqt(Q,dA); dA=1./dA;

U = 1 + 0*X;   %% Initial conditions
V = 1 + 0*X;
T = 0 + 0*X;
end

