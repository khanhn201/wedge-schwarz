function [Fx_total, Fy_total] = compute_drag_lift(U,V,P,nu,Xr,Rx,Dh,w, unxa_v, unya_v)
Fx_total = 0;
Fy_total = 0;
[dUdx,dUdy] = grad(U, Rx, Dh);
[dVdx,dVdy] = grad(V, Rx, Dh);
tau_xx = -P + 2*nu.*dUdx;
tau_xy =     nu.*(dUdy + dVdx);
tau_yy = -P + 2*nu.*dVdy;

Nelx = 16;

sidx = 1; % bottom wall
for r = 1:size(U, 1)
  wr = w(r);
  for e = 1:Nelx
    nx =  -unxa_v(r,e,sidx);
    ny =  -unya_v(r,e,sidx);
    % traction components at (r,e,sidx)
    T_x = tau_xx(r,e,sidx)*nx + tau_xy(r,e,sidx)*ny;
    T_y = tau_xy(r,e,sidx)*nx + tau_yy(r,e,sidx)*ny;
    Fx_total = Fx_total + T_x;
    Fy_total = Fy_total + T_y;
  end
end
