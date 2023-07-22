function [ exact_alpha ] = pickK_element_once( A, index )
%GAUSSIAN Summary of this function goes here
%   Detailed explanation goes here

% A matrix dimension
[m,n] = size(A);
kv = size(index,1);

% Compute exact alpha_kv
CtA=[];
for Cti=1:kv
   CtA=[-1*ones(2^(Cti-1),1),CtA;1*ones(2^(Cti-1),1),CtA];     
end
v_target = zeros(1,2^kv/2);
for v_index=1:2^(kv-1)
%     % CVX Solver
%     cvx_quiet true
%     cvx_begin
%         variable g0(n);
%         expression v0(m);
%         expression v_target0;
%         v0 = A*g0;
%         for z = 1:kv
%             v_target0 = v_target0 + CtA(v_index,z)*g0(index(z,1),1);
%         end
%         maximize( v_target0 );
%     subject to 
%         norm(g0,1) <= 1
%         0 <= v0 <= 0
%         0 <= (CtA(v_index,:)').*g0(index(:,1),1)
%     cvx_end
%     if(cvx_optval == Inf)
%         v_target0 = inf;
%     end
    
    indexS = CtA(v_index,:)'.*index;
    % linear programming
    f = zeros(2*n,1);
    f(abs(index(:,1)),1) = sign(indexS(:,1));
    b = [zeros(2*n+1,1);1];
    B = [A, zeros(m,n)];
    Aprime = [eye(n), -eye(n);-eye(n), -eye(n);-f';zeros(1,n),ones(1,n)];
    c = zeros(m,1);    
    [x,~,exitflag] = linprog(-f,Aprime,b,B,c); 
    if exitflag == -3 
        v = inf;
    end
    v = f'*x;
    v_target(1,v_index) = v;
end
exact_alpha = max(v_target);
end


  

