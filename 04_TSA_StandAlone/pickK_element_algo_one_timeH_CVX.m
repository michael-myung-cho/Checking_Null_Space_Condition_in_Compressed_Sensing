function [ a_portion ] = pickK_element_algo_one_timeH( H, sorted_alphaL_list, index, kv, L )
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
    v_target_exact = 0;
    cvx_begin
        variable g_exact(m);
        expression v_exact(n);
        expression v_target_exact;   
        v_exact = H*g_exact;
        vC_exact = v_exact;
        for z = 1:size(index,1)
           v_target_exact = v_target_exact + sign(index(z,1))*g_exact(abs(index(z,1)),1);
           vC_exact(abs(index(z,1)),1) = 0;
        end
        
        maximize( v_target_exact );
        
    subject to 
        norm(vC_exact,1) <= 1
        0 <= sign(index(:,1)).*v_exact(abs(index(:,1)),1)
    cvx_end
    if(cvx_optval == Inf)
        v_target_exact = inf;
    end
    a_portion = v_target_exact/(v_target_exact+1);    

end
end


  

