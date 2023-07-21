function [alpha,X,Y,Z]=BoundCVX(P,kv)
% Compute relaxation bound on alpha
[n,nk]=size(P);
disp('-------------------------------------------------------------------------------------');
disp('Solve SDP relaxation using CVX');
if kv>1
    nk=size(P,2);
    cvx_begin
    variable Q1(nk,nk) symmetric;
    variable Q2(n,nk);
    variable Q3(n,n) symmetric;
    variable t(n);
    variable r(n);
    maximize(trace(Q2'*P));
    subject to
        [Q1,Q2';Q2,Q3]==semidefinite(n+nk,n+nk)
        sum(sum(abs(P*Q1*P')))<=1
        max(abs(Q3))'<=t;
        abs(Q3)*ones(n,1)<=kv*t;
        sum(t)<=kv;
        max(t)<=1;
        max(abs(Q2*P'))'<=r;
        abs(P*Q2')*ones(n,1)<=kv*r;
        sum(r)<=1;
    cvx_end
    X=P*Q1*P';Y=Q3;Z=P*Q2';
    alpha=trace(Q2'*P);
else % When k=1, the extra columnwise constraints are redundant
    nn=n*(n+1)/2;
    nkn=(n+nk)*(n+nk+1)/2;
    clear blk;blk{1,1}='s';blk{1,2}=(n+nk);blk{1,3}=[2*ones(1,2*nn),2*ones(1,n^2)]; % Low rank SDP constraints
    At{1,1}=sparse(nkn,3);At{1,2}=[];At{1,3}=[];
    [Vb,Db]=eig([0 1; 1 0]);em=eye(n);
    for j=1:n;for i=1:j
            At{1,2}=[At{1,2},[P(i,:)',P(j,:)';zeros(n,2)]*Vb];
            if (i==j);At{1,3}=[At{1,3};[-1/2;1/2]];
            else At{1,3}=[At{1,3};[-1;1]];
            end; end
    end
    for j=1:n; for i=1:j
            At{1,2}=[At{1,2},[zeros(nk,2);em(:,i),em(:,j)]*Vb];
            if (i==j); At{1,3}=[At{1,3};[-1/2;1/2]];
            else; At{1,3}=[At{1,3};[-1;1]];
            end; end
    end
    for j=1:n; for i=1:n
            At{1,2}=[At{1,2},[P(i,:)',zeros(nk,1);zeros(n,1),em(:,j)]*Vb];
            At{1,3}=[At{1,3};[-1;1]*0.5];
        end
    end
    C{1,1}=[zeros(nk,nk),-0.5*P';-0.5*P,zeros(n,n)];
    blk{2,1}='l';blk{2,2}=4*nn+2*n^2; % Block linear constraints (switch U4 and U1)
    AA=[-sparse(ones(nn,1)),sparse(nn,2),speye(nn),sparse(nn,nn+n^2)];
    AA=[AA;-sparse(ones(nn,1)),sparse(nn,2),-speye(nn),sparse(nn,nn+n^2)];
    AA=[AA;sparse(nn,1),-sparse(ones(nn,1)),sparse(nn,1),sparse(nn,nn),speye(nn),sparse(nn,n^2)];
    AA=[AA;sparse(nn,1),-sparse(ones(nn,1)),sparse(nn,1),sparse(nn,nn),-speye(nn),sparse(nn,n^2)];
    AA=[AA;sparse(n^2,2),-sparse(ones(n^2,1)),sparse(n^2,2*nn),-speye(n^2)];
    AA=[AA;sparse(n^2,2),-sparse(ones(n^2,1)),sparse(n^2,2*nn),speye(n^2)];
    At{2,1}=AA; %clear AA;
    C{2,1}=sparse(4*nn+2*n^2,1);
    b=-[1;1;1;zeros(2*nn+n^2,1)]; % Objective
    OPTIONS.rmdepconstr = 1;
    %OPTIONS.gaptol=1e-2; %  Low target precision
    tic;[obj,M,y,Z]=sqlp(blk,At,C,b,OPTIONS);toc % Solve
    alpha=-b'*y;
    disp(['Result: ',num2str(-b'*y)])
    X=P*(M{1}(1:nk,1:nk))*P';Y=M{1}(nk+1:nk+n,nk+1:nk+n);Z=P*(M{1}(1:nk,nk+1:nk+n));
end
disp('-------------------------------------------------------------------------------------');