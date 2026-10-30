function [gap,nBad]=ls_audit(x,L,r,tol)
N=size(x,1);gap=inf;nBad=0;
for i=1:N-1
    dp=abs(x(i+1:end,:)-x(i,:));dp=min(dp,L-dp);
    gg=hypot(dp(:,1),dp(:,2))-2*r;
    gap=min(gap,min(gg));nBad=nBad+sum(gg < -tol);
end
end
