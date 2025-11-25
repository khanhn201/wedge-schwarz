function [xs, ys] = integrate_streamline(x0, y0, dt, nsteps, X, Y, U, V, P, T, z)
  xs = zeros(1, nsteps);
  ys = zeros(1, nsteps);
  xs(1) = x0;
  ys(1) = y0;

  for k = 2:nsteps
      x = xs(k-1);
      y = ys(k-1);

      % Get local velocity from your interpolate function
      [Uinterp, Vinterp, ~, ~] = interpolate([x], [y], X, Y, U, V, P, T, z);
      u1 = Uinterp;
      v1 = Vinterp;

      [Uinterp, Vinterp, ~, ~] = interpolate([x + 0.5*dt*u1], [y + 0.5*dt*v1], X, Y, U, V, P, T, z);
      u2 = Uinterp;
      v2 = Vinterp;

      [Uinterp, Vinterp, ~, ~] = interpolate([x + 0.5*dt*u2], [y + 0.5*dt*v2], X, Y, U, V, P, T, z);
      u3 = Uinterp;
      v3 = Vinterp;

      [Uinterp, Vinterp, ~, ~] = interpolate([x + dt*u3], [y + dt*v3], X, Y, U, V, P, T, z);
      u4 = Uinterp;
      v4 = Vinterp;

      % RK4 update
      xs(k) = x + dt/6 * (u1 + 2*u2 + 2*u3 + u4);
      ys(k) = y + dt/6 * (v1 + 2*v2 + 2*v3 + v4);
  end

