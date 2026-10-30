function H=ls_heap_create(cap)
H.t=zeros(cap,1);H.i=zeros(cap,1,'uint32');H.j=H.i;
H.vi=H.i;H.vj=H.i;H.nx=zeros(cap,1);H.ny=zeros(cap,1);
H.count=0;
end
