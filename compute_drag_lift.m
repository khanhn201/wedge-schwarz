function [Fx_total, Fy_total] = compute_drag_lift(U,V,P,nu,Xr,Rx,Dh,w)
sidx = 1; % bottom wall
Fx_total = 0;
Fy_total = 0;
[dUdx,dUdy] = grad(U, Rx, Dh);
[dVdx,dVdy] = grad(V, Rx, Dh);
tau_xx = -P + 2*nu.*dUdx;
tau_xy =     nu.*(dUdy + dVdx);
tau_yy = -P + 2*nu.*dVdy;

for r = 1:size(U, 1)
  wr = w(r);
  for e = 1:size(U, 2)
    tx = Xr(r,e,sidx,1,1); % dx/dr
    ty = Xr(r,e,sidx,2,1); % dy/dr
    tnorm = sqrt(tx^2 + ty^2);
    if tnorm == 0; continue; end
    nx =  ty/tnorm; ny = -tx/tnorm;   % outward unit normal candidate
    % traction components at (r,e,sidx)
    T_x = tau_xx(r,e,sidx)*nx + tau_xy(r,e,sidx)*ny;
    T_y = tau_xy(r,e,sidx)*nx + tau_yy(r,e,sidx)*ny;
    dS = tnorm * wr; % physical line element
    Fx_total = Fx_total + T_x * dS;
    Fy_total = Fy_total + T_y * dS;
  end
end
