function [permutated_A, sorted_alphaL_list, sorted_alphaR_list, sorted_alphaL_list_for_layerL] = Matrix_Permutation(A, kv, L) 

[m,n] = size(A);

% pick-l-step algo. and find remaining step algo. 
R = rem(kv,L);

% calculate alpha_1 from pick1_elementl algorithm
alpha1_list = pickL_element_algo_Lin(A,1);

% from the result of the pick1_element_algorithm, rearrange A matrix and obtain sorted alpha1 list 
A_alpha1 = zeros(n,m+1);
A_alpha1(:,m+1) = alpha1_list(:,2);
A_alpha1(:,1:m) = A';
permutated_A_alpha1 = -sortrows(-A_alpha1, m+1);
permutated_A = permutated_A_alpha1(:,1:m)';
alpha1_list(:,2) = permutated_A_alpha1(:,m+1);
sorted_alpha1_list = [ alpha1_list; -alpha1_list(:,1), alpha1_list(:,2)];
sorted_alpha1_list = -sortrows(-sorted_alpha1_list,2);

% from the permuted H matrix, obtain alpha_L list considered remaining part
% alphaL_list = pickL_element_algo(permutated_A, L)
% % Chain operation
% % m2 = size(alphaL_list,1);
% % alphaL_list(:,end+1) = 0;
% % for i=1:m2
% %     ind = 1;
% %     while(m2 - ind)
% %         if abs(alphaL_list(i,end-2)) < abs(alphaL_list(ind,1))
% %             alphaL_list(i,end) = alphaL_list(i,end-1) + alphaL_list(ind,end-1); 
% %             break;
% %         end
% %         ind = ind + 1;
% %     end
% % end
% sorted_alphaL_list = [ alphaL_list; -1*alphaL_list(:,1:end-2), alphaL_list(:,end-1:end)];
% sorted_alphaL_list = -sortrows(-sorted_alphaL_list, L+1);
% %sorted_alphaL_list(:,end) = [];
% sorted_alphaL_list_for_layerL = -sortrows(-alphaL_list, L+1);
% %sorted_alphaL_list_for_layerL(:,end) = [];


if L == 1
    sorted_alphaL_list = sorted_alpha1_list;
    sorted_alphaL_list_for_layerL = -sortrows(-alpha1_list, L+1); % only for layer L (reduce sign pattern cases: a{1,2,3} = a{-1,-2,-3}
else
    alphaL_list = pickL_element_algo_Lin(permutated_A, L);
    sorted_alphaL_list = [ alphaL_list; -alphaL_list(:,1:end-1), alphaL_list(:,end)];
    sorted_alphaL_list = -sortrows(-sorted_alphaL_list, L+1);
    sorted_alphaL_list_for_layerL = -sortrows(-alphaL_list, L+1);
end

% obtain alpha_r(remaining) list
if R == 1
    sorted_alphaR_list = sorted_alpha1_list;
elseif R == 0
    sorted_alphaR_list = 0;
else
    alphaR_list = pickL_element_algo_Lin(permutated_A, R); 
    sorted_alphaR_list = [ alphaR_list; -alphaR_list(:,1:end-1), alphaR_list(:,end)];
    sorted_alphaR_list = -sortrows(-sorted_alphaR_list, R+1);
end






