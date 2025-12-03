##function [xs, ys] = integrate_streamline(x0, y0, dt, nsteps, X, Y, X_tip, Y_tip, U,U_tip, V,V_tip,z, z_tip)
##  xs = zeros(1, nsteps);
##  ys = zeros(1, nsteps);
##  xs(1) = x0;
##  ys(1) = y0;
##
##  for k = 2:nsteps
##      x = xs(k-1);
##      y = ys(k-1);
##
##      [u1,v1] = interpolate_glo(x,y,X,Y,X_tip,Y_tip,U,U_tip,V,V_tip,z,z_tip);
##
##      [u2,v2] = interpolate_glo(x+0.5*dt*u1,y+0.5*dt*v1,X,Y,X_tip,Y_tip,U,U_tip,V,V_tip,z,z_tip);
##
##      [u3,v3] = interpolate_glo(x + 0.5*dt*u2,y + 0.5*dt*v2,X,Y,X_tip,Y_tip,U,U_tip,V,V_tip,z,z_tip);
##
##      [u4,v4] = interpolate_glo(x + dt*u3,y + dt*v3,X,Y,X_tip,Y_tip,U,U_tip,V,V_tip,z,z_tip);
##
##        % RK4 update
##        xs(k) = x + dt/6 * (u1 + 2*u2 + 2*u3 + u4);
##        ys(k) = y + dt/6 * (v1 + 2*v2 + 2*v3 + v4);
##
##        % Stop if outside domain or velocity is NaN
##        if isnan(xs(k)) || isnan(ys(k))
##            xs = xs(1:k-1);
##            ys = ys(1:k-1);
##            break;
##        end
##
##  end


function [xs, ys] = integrate_streamline(x0, y0,dt_init,nsteps, X, Y, X_tip, Y_tip, U, U_tip, V, V_tip, z, z_tip,C)
  xs = zeros(1, nsteps);
  ys = zeros(1, nsteps);
  xs(1) = x0;
  ys(1) = y0;

  % Estimate smallest grid spacing
  dx = mean(diff(unique(X(:))));
  dy = mean(diff(unique(Y(:))));
  dmin = min(dx, dy);

  for k = 2:nsteps
      x = xs(k-1);
      y = ys(k-1);

      % Get local velocity
      [u1, v1] = interpolate_glo(x, y, X, Y, X_tip, Y_tip, U, U_tip, V, V_tip, z, z_tip);
      speed = sqrt(u1^2 + v1^2);
      if speed == 0 || isnan(speed)
          xs = xs(1:k-1);
          ys = ys(1:k-1);
          break;
      end

      % Adaptive time step based on local velocity
      dt =  min(dt_init, C * dmin / speed);

      % RK4 integration
      [u2, v2] = interpolate_glo(x + 0.5*dt*u1, y + 0.5*dt*v1, X, Y, X_tip, Y_tip, U, U_tip, V, V_tip, z, z_tip);
      [u3, v3] = interpolate_glo(x + 0.5*dt*u2, y + 0.5*dt*v2, X, Y, X_tip, Y_tip, U, U_tip, V, V_tip, z, z_tip);
      [u4, v4] = interpolate_glo(x + dt*u3, y + dt*v3, X, Y, X_tip, Y_tip, U, U_tip, V, V_tip, z, z_tip);

      xs(k) = x + dt/6 * (u1 + 2*u2 + 2*u3 + u4);
      ys(k) = y + dt/6 * (v1 + 2*v2 + 2*v3 + v4);

      if isnan(xs(k)) || isnan(ys(k))
          xs = xs(1:k-1);
          ys = ys(1:k-1);
          break;
      end
  end
end

