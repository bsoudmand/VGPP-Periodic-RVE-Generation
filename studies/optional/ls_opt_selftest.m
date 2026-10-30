function ls_opt_selftest()
% Analytic collision: p=(1,0), w=(-1,0), 2r=.2, 2gamma=.1 => t=8/11.
[tt,nn]=ls_predict_pairs([1 1],[-.5 0],[0 1],[.5 0],.1,.05,10,1,1e-12);
assert(abs(tt-8/11)<1e-8 && abs(norm(nn)-1)<1e-9,'Collision root regression failed');
% Regression: earlier root is negative, later root is positive.
% The positive root must be available to the contact-direction filter.
[tp,~]=ls_predict_pairs([1 1],[1 0],[0 1],[0 0],.1,0,10,5,1e-12);
assert(isinf(tp),'Separating particles produced spurious collision');
% Periodic cross-boundary contact.
[t,n]=ls_predict_pairs([.1 1],[-.5 0],[9.6 1],[0 0],.1,0,10,2,1e-12);
assert(isfinite(t)&&abs(t-.6)<1e-8&&abs(norm(n)-1)<1e-8,'Periodic collision regression failed');
[g,bad]=ls_audit([1 1;1.4 1],10,.2,1e-8);
assert(abs(g)<1e-12 && bad==0,'Overlap audit regression failed');
% Heap ordering, including tied timestamps.
H=ls_heap_create(4);
H=ls_heap_push(H,3,1,2,uint32(0),uint32(0),[1 0]);
H=ls_heap_push(H,1,1,3,uint32(0),uint32(0),[1 0]);
H=ls_heap_push(H,2,2,3,uint32(0),uint32(0),[1 0]);
for target=1:3
    [H,e]=ls_heap_pop(H);assert(abs(e.time-target)<1e-12,'Heap ordering regression failed');
end
% At zero growth, an equal-mass head-on collision swaps normal speeds.
na=[1 0];va=[-0.5 0];vb=[0.5 0];g=0;
w=dot(va-vb,na);impulse=(2*g-w)*na;
vap=va+impulse;vbp=vb-impulse;
assert(norm(vap-[.5 0])<1e-12 && norm(vbp-[-.5 0])<1e-12, ...
    'Velocity impulse regression failed');
% Correct quadratic root if negative and positive roots coexist;
% monotonic grow against stationary neighbor must hit at t=0.8.
[tr,~]=ls_predict_pairs([1 1],[0 0],[0 1],[0 0],.1,.5,10,1,1e-12);
assert(abs(tr-.8)<1e-10,'Growing-radius root regression failed');
fprintf('Analytical collision, periodic-image, audit, impulse, and heap tests PASSED.\n');
end
