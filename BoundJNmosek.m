function alpha=BoundJNmosek(A,kv)
% Compute relaxation bound on alpha from Juditsky and Nemirovsky
disp('-------------------------------------------------------------------------------------');
disp('Solve Null Space Prop relaxation using J&N LP bound & MOSEK');
addpath('./JNcode');
m=size(A,1);
n=size(A,2);
% Solve for bound using MOSEK
rss.cs=InitiateS;
rss.m=size(A,1);
rss.n=size(A,2);
solver='M';
cm='C';
R=0.5*n/kv+sqrt(n);
rss.cs=GetBoundS(A,kv,'silent',rss.cs);
rss.fMP=1;
rss.fals=1;
alpha=rss.cs.gamma(kv);
disp(['Result: ',num2str(alpha)])
disp('-------------------------------------------------------------------------------------');