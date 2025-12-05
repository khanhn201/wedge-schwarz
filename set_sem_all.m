function [U,V,T,z,w,Dh,X,Y,Grr,Grs,Gss,Bl,Xr,Rx,Jac,Q,glo_num,Mu,Mv,Mp,Mt,ifnull,unxa_v,unya_v,BC_all,dA,interpdata_top]...
             = set_sem_all(N);

hdr;    % 2-D SEM multi-element

Nelx = 4;  Nely = 4; E = Nelx*Nely;
% Nelx = 1;  Nely = 1; E = Nelx*Nely;
N1=N+1;


%% Circular Geometry
x0 =  -1.5;  x1=-x0;  Lx=x1-x0;   % Domain coordinates and size
y0 =   1.0;  y1=1.5;  Ly=y1-y0;

%% Base Geometry
x0 =   0.0;  x1=2.;  Lx=x1-x0;   % Domain coordinates and size
y0 =  -0.5;  y1=1.5;  Ly=y1-y0;

zc = zwuni(Nelx); xc = x0 + Lx*(zc+1)/2;
zc = zwuni(Nely); yc = y0 + Ly*(zc+1)/2;

%% Problem parameters, as function of N

N1=N+1;
E=Nelx*Nely;

[z,w]=zwgll(N);                % Set basic operators
Dh=deriv_mat(z);

[R,S]=ndgrid(z,z);             % Build SEM mesh
X=zeros(N1,E,N1); Y=X;
alpha = (90-28.5/2) * pi/180;

e=0;
h = 0.35;
for ey=1:Nely; for ex=1:Nelx; e=e+1;
    function rb = bottom(r)
        slope = 0.075;
        if r <= 0.0
            rb = [sin(r*pi/2)*h*0.25/2, h*0.25*tan(alpha)/2 + slope*(sin(r*pi/2)+1)];
        else
            rb = [sin(r*pi/2)*h*0.25/2, h*0.25*tan(alpha)/2 + slope - slope*sin(r*pi/2)];
        end
        % rb = [r/2.0, 0.0];
    end
    function rb = top(r)
        rb = [1.0*(sin(r*pi/2))*0.25, 0.25*tan(alpha)];
        % rb = [1.0*(r+1.0)/2.0, 1.0];
    end
    function rb = blend(r, s)
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
##
##     node_id = 0;
##     for j = 1:N1; for i = 1:N1;
##       node_id = node_id +1;
##       text(Xe(i,j), Ye(i,j),num2str(node_id), 'fontsize',14);
##     end;end
##     xc = mean(Xe(:));
##     yc = mean(Ye(:));
##
##     % add element number
##     text(xc, yc, num2str(e), ...
##          'HorizontalAlignment','center', ...
##          'VerticalAlignment','middle', ...
##          'FontWeight','bold', ...
##          'Color','r');
##
## end
## pause;

 % [X,Y]=morph_circ(X,Y);         % Morph mesh

[Grr,Grs,Gss,Bl,Xr,Rx,Jac]=geom_elem(X,Y,Dh,w); % Terms for "A"
vol = sum(sum(sum(Bl)))
[Q,glo_num]=set_tp_semq(Nelx,Nely,N);


BC_all = [ 'D' 'D' 'D' 'D' ;     %% U
           'D' 'D' 'D' 'D' ;     %% V
           'N' 'N' 'N' 'N' ;     %% P
           'N' 'N' 'D' 'D' ];    %% T

[Mu,Q,glo_num]=set_mask(BC_all(1,:),Nelx,Nely,Q,glo_num);
[Mv,Q,glo_num]=set_mask(BC_all(2,:),Nelx,Nely,Q,glo_num);
[Mp,Q,glo_num]=set_mask(BC_all(3,:),Nelx,Nely,Q,glo_num); ifnull=1;
[Mt,Q,glo_num]=set_mask(BC_all(4,:),Nelx,Nely,Q,glo_num);

[unxa_v,unya_v] = set_unxy(Mu,X,Y,Xr);

dA=diag_sem(Grr,Grs,Gss,Dh); dA=qqt(Q,dA); dA=1./dA;

U = 0 - 0*X;   %% Initial conditions
V = 0 + 0*X;
T = 0 + 0*X;


% Store overlapping interpolation points
printf('**** Interpolation points at the top part *****\n')
interpdata_top = [];
n = 0;
for i = 1:Nelx
    for j = 1:size(X,1)
        n = n+1;
##        printf("Elem = %d, r = %f, s = %f, X = %f , Y = %f \n", i, j,1, X(j,i,1), Y(j,i,1));
        interpdata_top = [interpdata_top ; i, j,1, X(j,i,1), Y(j,i,1)];
    end
end
interpdata_top
##plot(X(:,1:Nelx,1), Y(:,1:Nelx,1), 'g')

##pause;

end

