function H=ls_heap_push(H,t,i,j,vi,vj,n)
if ~isfinite(t),return;end
if H.count>=numel(H.t)
    % Geometric array growth amortizes allocation cost.
    nNew=max(256,2*numel(H.t));
    H.t(nNew)=0;H.i(nNew)=uint32(0);H.j(nNew)=uint32(0);
    H.vi(nNew)=uint32(0);H.vj(nNew)=uint32(0);
    H.nx(nNew)=0;H.ny(nNew)=0;
end
k=H.count+1;H.count=k;
while k>1
    parent=floor(k/2);
    if H.t(parent)<=t,break;end
    H=ls_heap_copy(H,k,parent);k=parent;
end
H.t(k)=t;H.i(k)=uint32(i);H.j(k)=uint32(j);
H.vi(k)=vi;H.vj(k)=vj;H.nx(k)=n(1);H.ny(k)=n(2);
end
