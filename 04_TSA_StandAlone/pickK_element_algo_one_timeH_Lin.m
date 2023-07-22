function [ a_portion ] = pickK_element_algo_one_timeH_Lin( H, sorted_alphaL_list, index, kv, L )
%GAUSSIAN Summary of this function goes here
%   Detailed explanation goes here

% H matrix dimension: n*m 
[n,m] = size(H);
a_portion = 0;
if kv == L % index(1,1) small number, index(1,2) large number
    a_portion = find_portion(sorted_alphaL_list, index ,L);
else
    % ADMM Solver
%    a_portion = ADMM(A, index);
    
    % CVX Solver
%     v_target_exact = 0;
%     cvx_begin
%         variable g_exact(m);
%         expression v_exact(n);
%         expression v_target_exact;   
%         v_exact = H*g_exact;
%         vC_exact = v_exact;
%         for z = 1:size(index,1)
%            v_target_exact = v_target_exact + sign(index(z,1))*g_exact(abs(index(z,1)),1);
%            vC_exact(abs(index(z,1)),1) = 0;
%         end
%         
%         maximize( v_target_exact );
%         
%     subject to 
%         norm(vC_exact,1) <= 1
%         0 <= sign(index(:,1)).*v_exact(abs(index(:,1)),1)
%     cvx_end
%     if(cvx_optval == Inf)
%         v_target_exact = inf;
%     end
%     a_portion = v_target_exact/(v_target_exact+1);    

    % linear programming
    v = 0;
    f = zeros(n,1);
    f(abs(index(:,1)),1) = sign(index(:,1));
    b = [zeros(n,1);1;zeros(n+1,1)];
    I = eye(n);
    I(abs(index(:,1)),abs(index(:,1))) = 0;
    O = ones(1,n);
    O(1,abs(index(:,1))) = 0;
    IH = I*H; % n x m
    fTH = f'*H; % 1 x m
    fTHz = [fTH, zeros(1,n)];
    Aprime = [IH, -I;zeros(1,m), O;-IH, -I;-fTHz];
    [y,~,exitflag] = linprog(-fTHz,Aprime,b);
    if exitflag == -3
        v = inf;
    end
    v = fTHz*y;
    
    a_portion = v/(v+1);
end
end


  

