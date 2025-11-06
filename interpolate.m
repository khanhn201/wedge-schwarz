function [Uinterp,Vinterp,Tinterp] = interpolate(x,y,X,Y,U,V,T,z)
N1 = size(X,1);
E = size(X,2);

Uinterp = zeros(1, length(x));
Vinterp = zeros(1, length(x));
Tinterp = zeros(1, length(x));
for k = 1:length(x)
for e = 1:E
   % Corner coordinates of element e
   X11 = X(1,e,1);   Y11 = Y(1,e,1);
   X21 = X(N1,e,1);  Y21 = Y(N1,e,1);
   X12 = X(1,e,N1);  Y12 = Y(1,e,N1);
   X22 = X(N1,e,N1); Y22 = Y(N1,e,N1);

   % Coordinates of the element’s corners
   x_elem = [X11 X21 X22 X12];
   y_elem = [Y11 Y21 Y22 Y12];

   % Check if (x, y) is inside this element
   in = inpolygon(x(k), y(k), x_elem, y_elem);

   if in
     x_elem = [X11 X21 X22 X12];
     y_elem = [Y11 Y21 Y22 Y12];
     [r, s] = xy2rs(x(k), y(k), x_elem, y_elem);

      JR = interp_mat([r],z);
      JS = interp_mat([s],z);

      Uinterp(k) = tensor3(JS,1,JR,U(:,e,:));
      Vinterp(k) = tensor3(JS,1,JR,V(:,e,:));
      Tinterp(k) = tensor3(JS,1,JR,T(:,e,:));

##       fprintf('Point (%.3f, %.3f) is inside element %d | Uintp = %f, Vintp = %f, Tinterp = %f.\n', x(k), y(k), e, Uinterp(k),Vinterp(k),Tinterp(k));
       % You can now perform interpolation here...
       break
   end
end
end
end

function [r, s] = xy2rs(x, y, x_elem, y_elem)
% x_elem, y_elem = [X11 X21 X22 X12]
% initial guess
r = 0; s = 0;

for iter = 1:10
    % shape functions
    N = 0.25 * [(1 - r)*(1 - s);
                (1 + r)*(1 - s);
                (1 + r)*(1 + s);
                (1 - r)*(1 + s)];

    % derivatives wrt r and s
    dNdr = 0.25 * [-(1 - s);  (1 - s);  (1 + s); -(1 + s)];
    dNds = 0.25 * [-(1 - r); -(1 + r);  (1 + r);  (1 - r)];

    % compute current (x,y)
    x_r = sum(dNdr .* x_elem(:));
    y_r = sum(dNdr .* y_elem(:));
    x_s = sum(dNds .* x_elem(:));
    y_s = sum(dNds .* y_elem(:));

    x_cur = sum(N .* x_elem(:));
    y_cur = sum(N .* y_elem(:));

    % residual
    F = [x - x_cur; y - y_cur];

    % Jacobian matrix
    J = [x_r, x_s; y_r, y_s];

    % update (r,s)
    delta = J \ F;
    r = r + delta(1);
    s = s + delta(2);

    if norm(delta) < 1e-10
##      printf("norm(delta) = %f, not converged.\n", norm(delta))
##         resx = x - (0.25*(1-r)*(1-s)*x_elem(1) + 0.25*(1-r)*(s+1)*x_elem(4)+ 0.25*(r+1)*(1-s)*x_elem(2) + 0.25*(r+1)*(s+1)*x_elem(3) )
##         resy = y - (0.25*(1-r)*(1-s)*y_elem(1) + 0.25*(1-r)*(s+1)*y_elem(4)+ 0.25*(r+1)*(1-s)*y_elem(2) + 0.25*(r+1)*(s+1)*y_elem(3)  )
      return;
    end
end
##printf("norm(delta) = %f, not converged.\n", norm(delta))
end

