function [tau,nOut]=ls_predict_pairs(xa,va,xb,vb,r,gamma,L,horizon,epsT)
% Vectorized pair prediction for one particle a against many b, considering
% periodic images. Gives earliest approaching contact within horizon.
M=size(xb,1);tau=inf(M,1);nOut=zeros(M,2);
if M==0,return;end
w=va-vb; p0=xa-xb; R=2*r;G=2*gamma;
A=sum(w.^2,2)-G^2;
for sx=-1:1
    for sy=-1:1
        p=p0+L*[sx sy];
        B=2*(sum(p.*w,2)-R*G);
        C=sum(p.^2,2)-R^2;
        % Solve |p + w*t| = R + G*t using BOTH quadratic roots.
        % A smaller *negative* root must not hide a positive candidate.
        % Candidate is accepted only for an approaching contact.
        discr=B.^2-4*A.*C;
        cand=inf(M,1);
        lin=abs(A)<1e-13 & abs(B)>1e-13;
        tl=-C(lin)./B(lin);
        goodL=tl>=-epsT;
        idsL=find(lin);
        cand(idsL(goodL))=max(0,tl(goodL));
        quad=abs(A)>=1e-13 & discr>=0;
        idsQ=find(quad);
        if ~isempty(idsQ)
            sq=sqrt(max(0,discr(idsQ)));
            aq=A(idsQ);bq=B(idsQ);cq=C(idsQ);
            sg=ones(size(bq));sg(bq<0)=-1;
            q=-.5*(bq+sg.*sq);
            t1=q./aq;
            t2=inf(size(t1));
            stable=abs(q)>1e-30;
            t2(stable)=cq(stable)./q(stable);
            t2(~stable)=(-bq(~stable)+sq(~stable))./(2*aq(~stable));
            t1(t1 < -epsT)=inf;
            t2(t2 < -epsT)=inf;
            cand(idsQ)=max(0,min(t1,t2));
        end
        % Contact at current time with closing surface-normal speed.
        immediate=abs(C)<1e-10 & B<-1e-13;
        cand(immediate)=0;
        ok=isfinite(cand) & cand>=0 & cand<=horizon+1e-12 & cand<tau;
        ids=find(ok);
        if isempty(ids),continue;end
        pc=p(ids,:)+cand(ids).*w(ids,:);
        d=sqrt(sum(pc.^2,2));
        nn=pc./max(d,realmin);
        approach=sum(nn.*w(ids,:),2) < G-1e-11;
        ids=ids(approach);nn=nn(approach,:);
        tau(ids)=cand(ids);nOut(ids,:)=nn;
    end
end
end
