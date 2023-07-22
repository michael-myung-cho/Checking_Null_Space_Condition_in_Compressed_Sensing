function s=GetBoundS(A,sprs,mode,sold);
s=sold;
m=size(A,1);
n=size(A,2);
cls=m*n+n+1;
rws=2*n+1;
ctot=m*n+n*(n+1)+1;
rtot=n*rws;
if mode(1)=='?',
    go=input(sprintf('About to generate and solve %1dx%1d LP. To continue [y/n] > ',rtot,ctot),'s');
else
    disp(sprintf('Generating and solving %1dx%1d LP...',rtot,ctot));
    go='y';
end;
if go(1)~='y',
    gamma=inf;
    return;
end;
tstart=clock;
prb.a=[];
prb.blc=[];
prb.buc=[];
prb.blx=-inf*ones(m*n,1);
prb.bux=inf*ones(m*n,1);
base=0;
for i=1:n,
    a=A(:,i)';
    block=[];
    aft=ctot-base-cls;
    for irw=1:n,
        lsp=m*(irw-1);
        rsp=m*(n-irw)+irw-1;
        rw=[zeros(1,lsp),a,zeros(1,rsp),zeros(1,base),-1,zeros(1,n-irw),-1;...
            zeros(1,lsp),-a,zeros(1,rsp),zeros(1,base),-1,zeros(1,n-irw),-1];
        [ii,jj,vv]=find(rw);
        block=[block;sparse(ii,jj,vv,2,ctot)];
    end;
    rw=[zeros(1,m*n),zeros(1,base),ones(1,n),sprs,zeros(1,aft-1),-1];
    [ii,jj,vv]=find(rw);
    block=[block;sparse(ii,jj,vv,1,ctot)];
    blc=-inf*ones(rws,1);
    buc=zeros(rws,1);
    buc(2*i-1)=1;
    buc(2*i)=-1;
    buc(rws)=0;
    prb.a=[prb.a;block];
    prb.blc=[prb.blc;blc];
    prb.buc=[prb.buc;buc];
    prb.blx=[prb.blx;zeros(n,1);-inf];
    prb.bux=[prb.bux;inf*ones(n+1,1)];
    base=base+n+1;
end;
prb.blx=[prb.blx;0];
prb.bux=[prb.bux;inf];
prb.c=zeros(ctot,1);
prb.c(ctot)=1;
[rr,res]=mosekopt('minimize echo(0)',prb);
if rr>0,
    disp(sprintf('mosekopt cc: %1d',rr));
    if ~isfield(res,'rcode'),
        gamma=inf;
        Y=[];
        tused=etime(clock,tstart);
        disp(sprintf('Failure in computing alpha_s, CPU=%5.1f',tused));
        pause(0.1);
        return;
    end;
end;
if res.rcode>0,
    disp(res.rmsg);
end;
Y=zeros(m,n);
base=0;
for i=1:n,
    Y(:,i)=res.sol.itr.xx(base+1:base+m);
    base=base+m;
end;
Tmp=sort(abs(eye(n)-Y'*A));
s.gamma=[];
for i=1:size(Tmp,1);
    tmp=ones(1,i);
    gamma=max(tmp*Tmp(n-i+1:n,:));
    s.gamma=[s.gamma;gamma];
    if gamma<0.4999
        s.lp.als=max(s.lp.als,i);
    else
        if i>=sprs
            break;
        end
    end;
end;
s.lp.alsa=max(s.lp.als,s.lp.alsa);
s.lp.als=s.lp.alsa;
s.bst=SBest(s);
tused=etime(clock,tstart);
disp(sprintf('CPU=%5.1f',tused));
pause(0.01);