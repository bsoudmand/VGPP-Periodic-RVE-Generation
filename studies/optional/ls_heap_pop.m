function [H,e]=ls_heap_pop(H)
assert(H.count>0,'Cannot pop empty heap');
e.time=H.t(1);e.i=double(H.i(1));e.j=double(H.j(1));
e.vi=H.vi(1);e.vj=H.vj(1);e.normal=[H.nx(1) H.ny(1)];
last=H.count;H.count=last-1;
if last==1,return;end
v=H.t(last);ii=H.i(last);jj=H.j(last);a=H.vi(last);b=H.vj(last);
nx=H.nx(last);ny=H.ny(last);
k=1;
while 2*k<=H.count
    ch=2*k;
    if ch+1<=H.count && H.t(ch+1)<H.t(ch),ch=ch+1;end
    if H.t(ch)>=v,break;end
    H=ls_heap_copy(H,k,ch);k=ch;
end
H.t(k)=v;H.i(k)=ii;H.j(k)=jj;H.vi(k)=a;H.vj(k)=b;
H.nx(k)=nx;H.ny(k)=ny;
end
