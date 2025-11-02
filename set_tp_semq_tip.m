function [Q,glo_num]=set_tp_semq_tip(N);

hdr;

% Set up Q for tensor-product array of Nth-order elements in 2D

N1 = N+1;
E  = 1;
nL = E*N1*N1;

glo_num=zeros(N1,E,N1);

e=0;

n=0;

% Element1
 for j=1:N1; for i=1:N1;
    n=n+1;
    glo_num(i,1,j) = n;
 end; end;

## glo_num(:,2,1) = glo_num(:,1,N1);
## glo_num(1,3,:) = glo_num(N1,1,:);

## % Element2
## for j=2:N1; for i=1:N1;
##    n=n+1;
##    glo_num(i,2,j) = n;
## end; end;
## glo_num(:,3,N1) = glo_num(N1,2,:);
##
##% Element3
## for j=1:N1-1; for i=2:N1;
##    n=n+1;
##    glo_num(i,3,j) = n;
## end; end;

se_disp(glo_num,'glo_num')

Q=speye(nL);
k=0;
for j=1:N1;
for e=1:E;
for i=1:N1;
    ig = glo_num(i,e,j);
    k  = k+1;
    Q(k,k)=0; Q(k,ig)=1;
end;
end;
end;

Q=Q(:,1:n);
