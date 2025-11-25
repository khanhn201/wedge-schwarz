function [U,V,T,z,w,Dh,X,Y,Grr,Grs,Gss,Bl,Xr,Rx,Jac,Q,glo_num,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,BC_all,dA,interpdata_tip]...
             = set_sem_all_tipv03(N);

hdr;    % 2-D SEM multi-element

Nelx = 2;  Nely = 2; E = Nelx*Nely;
% Nelx = 1;  Nely = 1; E = Nelx*Nely;
N1=N+1;

%% Problem parameters, as function of N

N1=N+1;
##E=Nelx*Nely;

[z,w]=zwgll(N);                % Set basic operators
Dh=deriv_mat(z);

[R,S]=ndgrid(z,z);             % Build SEM mesh
X=zeros(N1,E,N1); Y=X;

% slopes from apex (constant on both sides)
alpha = (90-28.5/2) * pi/180;

e=0;

uf = @(x) -2*cos(pi/4*(x+1))+1;
h = 0.4;
for ey=1:Nely; for ex=1:Nelx; e=e+1;
    function rb = bottom(r)
        slope = 0.5;
        r = uf(r);
        rb = [(r+1.0)/2.0*h*0.25/2, ((r+1)/2)*h*0.25*tan(alpha)/2];

    end
    function rb = top(r)
        r = uf(r);
        rb = [ (r-1)/2*h*0.25/2, (0.02*(r+1)/2)*tan(alpha) + h*0.25*tan(alpha)/2];
    end
    function rb = blend(r, s)
        s = uf(s);
        rb = top(r)*(s+1.0)/2.0 + bottom(r)*(-s+1.0)/2.0;
    end

    r00 = blend((ex-1)/Nelx*2-1, (ey-1)/Nely*2-1);
    r10 = blend((ex)/Nelx*2-1, (ey-1)/Nely*2-1);
    r01 = blend((ex-1)/Nelx*2-1, (ey)/Nely*2-1);
    r11 = blend((ex)/Nelx*2-1, (ey)/Nely*2-1);

    xe00 = r00(1);  ye00 = r00(2);
    xe10 = r10(1);  ye10 = r10(2);
    xe01 = r01(1);  ye01 = r01(2);
    xe11 = r11(1);  ye11 = r11(2);

    X(:,e,:) = (xe00 + (xe10-xe00).*(R+1)/2) .* (-S+1)/2 ...
             + (xe01 + (xe11-xe01).*(R+1)/2) .* (1+S)/2;
    Y(:,e,:) = (ye00 + (ye10-ye00).*(R+1)/2) .* (-S+1)/2 ...
             + (ye01 + (ye11-ye01).*(R+1)/2) .* (1+S)/2;
end; end;

## figure; hold on; axis equal;
## for e = 1:E
##     % Extract patch for element e
##     Xe = squeeze(X(:,e,:));
##     Ye = squeeze(Y(:,e,:));
##     % Plot grid lines
##     plot(Xe, Ye, 'k-');           % lines along R
##     plot(Xe', Ye', 'k-');         % lines along S

##     node_id = 0;
##     for j = 1:N1; for i = 1:N1;
##       node_id = node_id +1;
##       text(Xe(i,j), Ye(i,j),num2str(node_id), 'fontsize',14);
##     end;end
##     xc = mean(Xe(:));
##     yc = mean(Ye(:));

##     % add element number
##     text(xc, yc, num2str(e), ...
##          'HorizontalAlignment','center', ...
##          'VerticalAlignment','middle', ...
##          'FontWeight','bold', ...
##          'Color','r');

## end
## pause;

 % [X,Y]=morph_circ(X,Y);         % Morph mesh

[Grr,Grs,Gss,Bl,Xr,Rx,Jac]=geom_elem(X,Y,Dh,w); % Terms for "A"
vol = sum(sum(sum(Bl)))
[Q,glo_num]=set_tp_semq(Nelx,Nely,N);


BC_all = [ 'D' 'D' 'D' 'D' ;     %% U
           'D' 'D' 'D' 'D' ;     %% V
           'N' 'D' 'N' 'D' ;     %% P
           'N' 'D' 'N' 'D' ];    %% T

[Mu,Q,glo_num]=set_mask(BC_all(1,:),Nelx,Nely,Q,glo_num);
[Mv,Q,glo_num]=set_mask(BC_all(2,:),Nelx,Nely,Q,glo_num);
[Mp,Q,glo_num]=set_mask(BC_all(3,:),Nelx,Nely,Q,glo_num); ifnull=1;
[Mt,Q,glo_num]=set_mask(BC_all(4,:),Nelx,Nely,Q,glo_num);

[unxa_v,unya_v] = set_unxy(Mu,X,Y,Xr);

dA=diag_sem(Grr,Grs,Gss,Dh); dA=qqt(Q,dA); dA=1./dA;

U = 0 + 0*X;   %% Initial conditions
V = 0 + 0*X;
T = 0 + 0*X;


% Store overlapping interpolation points
printf('**** Interpolation points at the tip part *****\n')
interpdata_tip = [];
n = 0;
for k = E-Nelx+1:E
   j = N1
   for i = 1:size(X,1)
       n = n+1;
       interpdata_tip = [interpdata_tip ; k, i,j, X(i,k,j), Y(i,k,j)];
   end
end
for k = Nelx:Nelx:E
   i = N1
   for j = 1:size(Y,1)
       n = n+1;
       interpdata_tip = [interpdata_tip ; k, i,j, X(i,k,j), Y(i,k,j)];
   end
end



interpdata_tip
##plot(X(N1,2,:), Y(N1,2,:), 'g')
##hold on;
##plot(X(:,3,N1), Y(:,3,N1), 'g')
##plot(X(N1,4,:), Y(N1,4,:), 'g')
##plot(X(:,4,N1), Y(:,4,N1), 'g')
##hold off;
##size(interpdata_tip)
##pause
end

