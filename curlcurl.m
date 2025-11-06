function [curlcurlX,curlcurlY,Omega]=curlcurl(U,V,Bl,Rx,Dh);
%%   Evaluate curl-curl term (to be extrapolated)

     % Omega     = tensor3(1,1,Dh,V)-tensor3(Dh,1,1,U);
     % curlcurlX =  Bl.*tensor3(Dh,1,1,Omega);
     % curlcurlY = -Bl.*tensor3(1,1,Dh,Omega);


R11 = Rx(:,:,:,1,1); R12 = Rx(:,:,:,1,2);
R21 = Rx(:,:,:,2,1); R22 = Rx(:,:,:,2,2);

Ur = tensor3(1,1,Dh,U);
Us = tensor3(Dh,1,1,U);
Vr = tensor3(1,1,Dh,V);
Vs = tensor3(Dh,1,1,V);
Ux = R11 .* Ur + R21 .* Us;
Uy = R12 .* Ur + R22 .* Us;
Vx = R11 .* Vr + R21 .* Vs;
Vy = R12 .* Vr + R22 .* Vs;

Omega = Vx - Uy;

Or = tensor3(1,1,Dh,Omega);
Os = tensor3(Dh,1,1,Omega);
Ox = R11 .* Or + R21 .* Os;
Oy = R12 .* Or + R22 .* Os;

curlcurlX =  Bl .* Oy;
curlcurlY = -Bl .* Ox;
