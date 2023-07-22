function [snew,F,tused]=GetAlpha1Bound(A,sold,ctrl)
%GetAlpha1Bound computes a lower bound for the supported by A sparsity
%
%   [s_new,Y]=GetAlpha1Bound(A,s_old,control);
%
%   Inputs:
%   A; (double array) matrix to be tested
%   printlevel: (0 or 1) controls the information displayed during the
%   execution. 
%   s_old: the structure of results. Should contain the fields:
%       s_mi: sparsity bound obtained using the mutual incoherence
%       s_bst: the best sparsity bound obtained
%       s_alpha: structure containing sparsity bounds 
%       Fields of s_alpha:
%           al1: an "exact" upper bound on alpha_1 (as in the J.-N.'s paper),
%           al1s: an upper bound on alpha_s obtained using the matrix eye(n)-Y'*A
%           computed in course of alpha_1 computation
%           s_alpha_1: corresponding sparsity bound
%           als: alpha_s bound
%           s_alpha_s: corresponding sparsity bound
%       control: structure with parameters of the procedure
%       Fields of control:
%           solver: either 'LP' (or 'lp', etc.) or 'MP' (or 'mp', or
%           'Mirror',etc.) -- the solver used. LP reformulation is soved
%           using MOSEK; saddle-point reformulation is solved using
%           Mirror-Prox algorithm
%           beta: [1e6] bound on the Euclidean norm of columns of the matrix certificate 
%           
%   Outputs: 
%   s_new: new structure (same as s_old) containing the updated sparsity bounds
%   Y: matrix certificate (matrix of the same dimension as A) as in the J.-N.'s paper
tstart=clock;
ni=nargin;

%if ~ni, error('At least one input argument is required!'); end
[m,n]=size(A);
printlevel=1;
half=0.499999;
cond=1e-8;
solver=0; % Using Mosek by default 
beta=1e6;
if ni>=3, 
    if isfield(ctrl,'printlevel'), printlevel=ctrl.printlevel; end
    if isfield(ctrl,'half'), half=ctrl.half; end
    if isfield(ctrl,'cond'), cond=ctrl.cond; end
    if isfield(ctrl,'beta'), beta=ctrl.beta; end

    if isfield(ctrl,'solver'), 
        if ~isempty(strfind(lower(ctrl.solver),'m')), solver=1; end% MP-solver
    end
end
%
if ni<2||isempty(sold)
    sold.s_alpha.al1=inf;
    sold.s_alpha.al1s=inf;
    sold.s_alpha.als=inf;
    sold.s_alpha.s_alpha_1=0;
    sold.s_alpha.s_alpha_s=0;
    sold.s_mi=0; 
end
snew=sold;
% compute mutual incoherehce
if isempty(snew.s_mi)|| snew.s_mi==inf 
    gmx=0;
    for i=1:n,
    y=A(:,i);
    y=y/norm(y)^2;
    z=A'*y;
    z(i)=z(i)-1;
    g=max(abs(z));
    g=g/(1+g);
    gmx=max(gmx,g);
    end;
    snew.s_mi=max(min(floor(half/gmx),m),snew.s_mi);
    if printlevel,
        fprintf('s.s_mi=%1d',snew.s_mi);
    	drawnow
    end;
end
% "normalize" sensing matrix
[U,D,V]=svd(A);
dm=max(diag(D));
I=find(diag(D)>cond*dm);
B=V(:,I)';
if solver,
    [alpha1,Vs,Valpha,Y]=GetAlpha1BoundMP(B,half,printlevel);
else
    [alpha1,Vs,Valpha,Y]=GetAlpha1BoundLP(B,beta,half,printlevel);
end

snew.s_alpha.al1=min(snew.s_alpha.al1,max(alpha1));
snew.s_alpha.al1s=min(snew.s_alpha.al1s,Valpha);
snew.s_alpha.s_alpha_1=max(snew.s_alpha.s_alpha_1, Vs);
U1=U(:,I); D1=diag(1./diag(D(I,I)));
F=U1*D1*Y;
snew.bst=max([snew.s_alpha.s_alpha_s,snew.s_alpha.s_alpha_1,snew.s_mi]);

tused=etime(clock,tstart);
if printlevel,
    fprintf('Certified sparsity: %d (Alpha_%d=%2.7f); ',Vs,Vs,Valpha);
    fprintf('CPU=%5.1f  \n',tused);
    drawnow
end
end
%
function [Ral,Vs,Valpha,Y]=GetAlpha1BoundLP(A,beta,half,prn)
%Ral: computed values of alpha_1
%Vs: the maximum sparsity that is verified using the solution found
%Valpha: the value of alpha corresponding to Vs verification
Vs=0;
Valpha=0;
[m,n]=size(A);
Ral=zeros(n,1);
gmx=0;
%
[ii,jj,vv]=find([[A';-A'],-ones(2*n,1),zeros(2*n,m);eye(m),zeros(m,1),-eye(m);-eye(m),zeros(m,1),-eye(m);zeros(1,m+1),ones(1,m)]);
prb.a=sparse(ii,jj,vv,2*n+2*m+1,2*m+1);
prb.blc=-inf*ones(2*n+2*m+1,1);
prb.buc=zeros(2*n+2*m+1,1);
prb.buc(end)=beta;
prb.blx=-inf*ones(2*m+1,1);
prb.bux=inf*ones(2*m+1,1);
prb.c=zeros(2*m+1,1);
prb.c(m+1)=1;

clear ii;
clear jj;
clear vv;
nvar=size(prb.a,2);
ncns=size(prb.a,1);
nnul=nnz(prb.a);
fprintf('Using MOSEK to solve  %1d %1dx%1d LPs with %1d nonzeros \n',n,ncns,nvar,nnul);
Y=zeros(m,n);
%
if prn,
    oped=0;
    fprintf(1,'%s',sprintf('  0%% done')); drawnow; 
end
for i=1:n,
    if (i>1),
        prb.buc(i-1)=0;
        prb.buc(n+i-1)=0;
    end;
    prb.buc(i)=1;
    prb.buc(i)=1;
    prb.buc(n+i)=-1;
    if prn>1,
        [rr,res]=mosekopt('minimize',prb);
    else
        [rr,res]=mosekopt('minimize echo(0)',prb);
    end
    if rr>0,
        if isfield(res,'rcode')&&(res.rcode>0), error(res.rmsg);
        else error(['mosekopt rr: ',num2str(rr)]);
        end;
    end
    Ral(i)=prb.c'*res.sol.itr.xx;
    Y(:,i)=res.sol.itr.xx(1:m);
    tmp=res.sol.itr.xx(1:m)'*A;
    if prn
        ped=floor(100*i/n);
        if ped>oped,
            fprintf(1,'%s',...
                sprintf('\b\b\b\b\b\b\b\b\b%3u%% done',ped));
            drawnow
            oped=ped;
        end
    end

end;
Vs=floor(half/max(Ral));
Tmp=sort(abs(eye(n)-Y'*A),1,'ascend');
for sc=1:n,
    tmp=sum(Tmp(n-sc+1:n,:),1);
    if max(tmp)<half,
        Vs=max(Vs,sc);
        if Vs==sc,
            Valpha=max(tmp);
        end;
    else
        break;
    end;
end;
if prn, fprintf(1,'%s',sprintf('\b\b\b\b\b\b\b\b\b')); drawnow; end
end
%
function [Ral,Vs,Valpha,F]=GetAlpha1BoundMP(B,half,prn)
%Vs: the maximum sparsity that is verified using the solution found
%Valpha: the value of alpha corresponding to Vs verification
Vs=0;
Valpha=0;
[no,ns]=size(B);
if prn
    fprintf('Solving %1dx%1d saddle point reformulation of %1dx%1d LP... \n',2*ns*ns,no*ns,ns*ns+1,no*ns);
    drawnow;
end
[sol,F]=MPAlpha1(B,half);
slb=sol.lb;
sub=sol.ub;
cpu=sol.cpu;

Y=F'*B;
sY=sort(abs(eye(ns)-Y),1);
Ral=sY(ns,:)';
Vs=0;
for sc=1:ns,
    tmp=sum(sY(ns-sc+1:ns,:),1);
    if max(tmp)<half,
        Vs=max(Vs,sc);
        if Vs==sc,
            Valpha=max(tmp);
        end;
    else
        break;
    end;
end;
end