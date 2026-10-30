function H=ls_heap_copy(H,d,s)
H.t(d)=H.t(s);H.i(d)=H.i(s);H.j(d)=H.j(s);
H.vi(d)=H.vi(s);H.vj(d)=H.vj(s);
H.nx(d)=H.nx(s);H.ny(d)=H.ny(s);
end
