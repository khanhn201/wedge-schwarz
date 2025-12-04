clear; close all; clc;

##load('n30.mat')

##xplt = {}; yplt = {};
##
##% point 1
##for i=1:4
## x0 = 0.0;   % starting x
## y0 = 0.03+i/5*0.03;   % starting y
##
##
## C = 10;
## dtstream = 1e9;
## % Integrate streamline using your function
## [xs, ys] = integrate_streamline(x0, y0,dtstream, 200, X, Y, X_tip, Y_tip, U,U_tip, V,V_tip,z, z_tip,C);
##
## % Plot the streamline on top
## plxs 'b-', 'LineWidth', 1.5); hold on;drawnow;
##
## xplt{i} = xs;
## yplt{i} = ys;
##end
##
##save data_sl_loop1.mat

##load data_sl_loop1.mat
##
##for i = 1:4
##  plot(xplt{i},yplt{i},'b--','LineWidth',1.25); hold on;drawnow;
##end
##
##% point 1
##for i=1:4
## x0 = 0.0;   % starting x
## y0 = 0.06+i/10*0.1;   % starting y
##
## C = 100;
## dtstream = 1e8;
## % Integrate streamline using your function
## [xs, ys] = integrate_streamline(x0, y0,dtstream, 100, X, Y, X_tip, Y_tip, U,U_tip, V,V_tip,z, z_tip,C);
##
## % Plot the streamline on top
## plot(xs, ys, 'r-', 'LineWidth', 1.5); hold on;drawnow;
## xplt{i+4} = xs;
## yplt{i+4} = ys;
##end
##
##save data_sl_loop2.mat

##load data_sl_loop2.mat
##for i = 1:8
##  plot(xplt{i},yplt{i},'b--','LineWidth',1.25); hold on;drawnow;
##end
##
##% point2
##for i=1:5
## x0 = 0.0;   % starting x
## y0 = 0.13+(i-1)/4*0.1;   % starting y
##
## dtstream = 1e4;
## C = 100;
## % Integrate streamline using your function
## [xs, ys] = integrate_streamline(x0, y0, dtstream, 200, X, Y, X_tip, Y_tip, U,U_tip, V,V_tip,z, z_tip,C);
##
## % Plot the streamline on top
## plot(xs, ys, 'r-', 'LineWidth', 1.5); hold on;drawnow;
##  xplt{i+8} = xs;
##  yplt{i+8} = ys;
##end
##
##save data_sl_loop3.mat

##load data_sl_loop3.mat
##for i = 1:13
##  plot(xplt{i},yplt{i},'b-','LineWidth',1.25); hold on;drawnow;
##end
##
##% point3
##for i=1:4
## x0 = 0.0;   % starting x
## y0 = 0.26+(i-1)/4*0.3;   % starting y
##
## dtstream = 1e4;
## C = 100;
## % Integrate streamline using your function
## [xs, ys] = integrate_streamline(x0, y0, dtstream, 300, X, Y, X_tip, Y_tip, U,U_tip, V,V_tip,z, z_tip,C);
##
## % Plot the streamline on top
## plot(xs, ys, 'r-', 'LineWidth', 1.5); hold on;drawnow;
## xplt{i+13} = xs;
## yplt{i+13} = ys;
##
##end
##
##save data_sl_loop4.mat

##load data_sl_loop4.mat
##for i = 1:17
##  plot(xplt{i},yplt{i},'b-','LineWidth',1.25); hold on;drawnow;
##end
##
##% point3
##for i=1:5
## x0 = 0.0;   % starting x
## y0 = 0.60+(i-1)/4*0.35;   % starting y
##
## dtstream = 1;
## C = 100;
## % Integrate streamline using your function
## [xs, ys] = integrate_streamline(x0, y0, dtstream, 500, X, Y, X_tip, Y_tip, U,U_tip, V,V_tip,z, z_tip,C);
##
## % Plot the streamline on top
## plot(xs, ys, 'r-', 'LineWidth', 1.5); hold on;drawnow;
## xplt{i+17} = xs;
## yplt{i+17} = ys;
##
##end
##
##save data_sl_loop5.mat

load data_sl_loop5.mat
for i = 1:21
  plot(xplt{i},yplt{i},'b-','LineWidth',1.25); hold on;drawnow;
end

alpha = (28.5/2) * pi/180;

x = [0, 1*tan(alpha), -1*tan(alpha), 0];
y = [0 , 0.984377 , 0.984377, 0];

plot(x,y,'k','LineWidth',1.2);


xlim([-0.5,0.5]);
ylim([0,1.2]);
axis equal;




