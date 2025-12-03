function [u1,v1] = interpolate_glo(x,y,X,Y,X_tip,Y_tip,U,U_tip,V,V_tip,z,z_tip)

N1 = size(X,1);
  % Corner coordinates of element e
      X11 = X_tip(1,1,1);   Y11 = Y_tip(1,1,1);
      X21 = X_tip(N1,2,1);  Y21 = Y_tip(N1,2,1);
      X12 = X_tip(1,3,N1);  Y12 = Y_tip(1,3,N1);
      X22 = X_tip(N1,4,N1); Y22 = Y_tip(N1,4,N1);

      % Coordinates of the element’s corners
      x_elem = [X11 X21 X22 X12];
      y_elem = [Y11 Y21 Y22 Y12];
      intip = inpolygon(x, y, x_elem, y_elem);


      if intip
        % Get local velocity from your interpolate function
        [u1, v1, ~, ~] = interpolate([x], [y], X_tip, Y_tip, U_tip, V_tip, 0*U_tip, 0*U_tip, z_tip);
      else
        [u1, v1, ~, ~] = interpolate([x], [y], X, Y, U, V, 0*U, 0*U, z);
      end

 end
