function [ exact_alphaK ] = exact_alphaK_algo( permutated_A, permutated_H, sorted_alphaL_list, best_candidate_node, kv, L )
%EXACT_ALPHAK_ALGO Summary of this function goes here

% find indexs of parent nodes from the tree structure in order
% to calculate the exact alphaK under given indexs
current_node = best_candidate_node;
current_layer = current_node.Property8;
index_for_exact_alpha = [];

% find indices
number_of_layers = ceil(current_layer/L);
for i_in = 1:number_of_layers
    index_for_exact_alpha = [current_node.Property7, index_for_exact_alpha];
    current_node = current_node.Prev;
end
% when L=2, k=5
if size(index_for_exact_alpha,2) > 1
    while index_for_exact_alpha(:,end) == index_for_exact_alpha(:,end-1)
        index_for_exact_alpha(:,end) = [];
    end
end
% pickk-element-algo
% index_for_exact_alpha: low to high absolute values order
if (current_layer <= L) && (current_layer < kv)
    exact_alphaK = best_candidate_node.Property1;
else
    [m,n] = size(permutated_A);
    if 2*m <= n
        exact_alphaK = pickK_element_algo_one_time_Lin(permutated_A, sorted_alphaL_list, index_for_exact_alpha(1,:)', kv, L );
    else
        exact_alphaK = pickK_element_algo_one_timeH_Lin(permutated_H, sorted_alphaL_list, index_for_exact_alpha(1,:)', kv, L );
    end
end

