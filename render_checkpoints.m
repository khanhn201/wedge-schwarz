hf = figure('Position',[100 100 1920 1080]);
for file_id=1:84
   fname = sprintf('ar05/c%06d.mat', file_id);
   load(fname);
   x_all = [];
   y_all = [];
   u_all = [];
   v_all = [];
   N_top = size(X,1);
   for e = 1:size(X,2)
       x = reshape(X(:,e,:), N_top, N_top);
       y = reshape(Y(:,e,:), N_top, N_top);
       u = reshape(U(:,e,:), N_top, N_top);
       v = reshape(V(:,e,:), N_top, N_top);
       x_all = [x_all; x(:)];
       y_all = [y_all; y(:)];
       u_all = [u_all; u(:)];
       v_all = [v_all; v(:)];
   end
   N_tip = size(X_tip,1);
   for e = 1:size(X_tip,2)
       x = reshape(X_tip(:,e,:), N_tip, N_tip);
       y = reshape(Y_tip(:,e,:), N_tip, N_tip);
       u = reshape(U_tip(:,e,:), N_tip, N_tip);
       v = reshape(V_tip(:,e,:), N_tip, N_tip);
       x_all = [x_all; x(:)];
       y_all = [y_all; y(:)];
       u_all = [u_all; u(:)];
       v_all = [v_all; v(:)];
   end
   mag = sqrt(u_all.^2 + v_all.^2);
   mag(mag == 0) = 1;  % avoid division by zero
   % scatter(x_all, y_all, 50, mag, 'filled'); hold on;
   % colormap(jet);
   % quiver(x_all, y_all, u_all, v_all, 'k'); hold off;
   hold off;
   
nx = 400; ny = 400;
xg = linspace(min(x_all), max(x_all), nx);
yg = linspace(min(y_all), max(y_all), ny);
[Xg, Yg] = meshgrid(xg, yg);

% Interpolate scattered data to grid
MagG = griddata(x_all, y_all, mag, Xg, Yg);

% Line contours
contourf(Xg, Yg, MagG, 40, 'LineColor', 'none');
colormap(jet);
colorbar;
axis equal tight;
hold on;
E1 = size(X,2);
for e = 1:E1
    Xe = squeeze(X(:,e,:));
    Ye = squeeze(Y(:,e,:));
    plot(Xe,  Ye,  'k-', 'LineWidth', 0.5);   % R-lines
    plot(Xe', Ye', 'k-', 'LineWidth', 0.5);   % S-lines
end

% ====== OVERLAY TIP GRID ======
E2 = size(X_tip,2);
for e = 1:E2
    Xe = squeeze(X_tip(:,e,:));
    Ye = squeeze(Y_tip(:,e,:));
    plot(Xe,  Ye,  'k-', 'LineWidth', 0.5);
    plot(Xe', Ye', 'k-', 'LineWidth', 0.5);
end

% Velocity vectors
% quiver(x_all, y_all, u_all, v_all, 'k');
hold off;

   fname = sprintf('ar05/c%06d.png', file_id);
   print (hf, fname);
end
