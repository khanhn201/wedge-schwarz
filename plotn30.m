load n30.mat

nn = 100;
x = zeros(nn);
y = linspace(0.98,0.005,nn);

u1plt = [];

for i = 1:nn
  [u1, v1] = interpolate_glo(x(i), y(i), X, Y, X_tip, Y_tip, U, U_tip, V, V_tip, z, z_tip);

  u1plt = [u1plt; u1];

endfor

figure;
plot(y, u1plt)

k = 0;umax = 1;
eddy_strength = [];
eddy_size = [];
y_n = 1;
for i = 2:nn;
  u0 = u1plt(i-1);
  u1 = u1plt(i);

  if u0*u1 <= 0.0; k= k+1
     eddy_strength(k) = umax;
     eddy_size(k) = y_n-y(i);
     umax = 0;
     y_n = y(i);
  end
  umax = max(umax,abs(u1));
end

figure;
semilogy(1:k, eddy_strength);


figure;
semilogy(1:k, eddy_size);
title("Eddy size")


figure;
plot(1:k-1,eddy_strength(1:k-1)./eddy_strength(2:k))

figure;
plot(1:k-1,eddy_size(1:k-1)./eddy_size(2:k))

