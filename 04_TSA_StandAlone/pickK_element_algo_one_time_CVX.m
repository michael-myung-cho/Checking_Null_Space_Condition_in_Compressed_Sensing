function [ a_portion ] = pickK_element_algo_one_time( A, sorted_alphaL_list, index, kv, L )
%GAUSSIAN Summary of this function goes here
%   Detailed explanation goes here

% H matrix dimension: n*m 
[m,n] = size(A);
a_portion = 0;
if kv == L % index(1,1) small number, index(1,2) large number
    a_portion = find_portion(sorted_alphaL_list, index ,L);
else
    % ADMM Solver
%    a_portion = ADMM(A, index);
    
    % CVX Solver
%     v_target_exact = 0;
%     cvx_begin
%         variable g_exact(n);
%         expression v_exact(m);
%         expression v_target_exact;   
%         v_exact = A*g_exact;
%         for z = 1:size(index,1)
%            v_target_exact = v_target_exact + sign(index(z,1))*g_exact(abs(index(z,1)),1);
%         end           
%         maximize( v_target_exact ); 
%     subject to 
%         norm(g_exact,1) <= 1
%         0 <= v_exact <= 0
%         0 <= sign(index(:,1)).*g_exact(abs(index(:,1)),1)
%     cvx_end
%     if(cvx_optval == Inf)
%         v_target_exact = inf;
%     end

    % linear programming
    f = zeros(2*n,1);
    f(abs(index(:,1)),1) = sign(index(:,1));
    b = [zeros(2*n+1,1);1];
    B = [A, zeros(m,n)];
    Aprime = [eye(n), -eye(n);-eye(n), -eye(n);-f';zeros(1,n),ones(1,n)];
    c = zeros(m,1);    
     [x,~,exitflag] = linprog(-f,Aprime,b,B,c); 
     if exitflag == -3 
         v = inf;
     end
     v = f'*x;
    a_portion = v;    

end
end


  

