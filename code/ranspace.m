function y=ranspace(X1,X2,N)
%ranspace RAK 28Sep15 (scalar inputs)
rr=cumsum(rand(1,N+2));
g=(X2-X1)/(rr(end)-rr(1));
y=g*(rr-rr(1))+X1;
y=y(2:end-1);
