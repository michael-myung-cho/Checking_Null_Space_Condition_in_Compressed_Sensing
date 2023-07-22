function [ node, flag_ok_to_make ] = make_one_best_candidate_node( target_node, node, sorted_alphaL_list, sorted_alphaR_list, sorted_alphaL_list_for_layerL, k, L, max_index, max_upper_bound_layer)
% MAKE_ONE_BEST_CANDIDATE_NODE is function for attaching
% one_best_candidate_node in the tree structure.
% (NOTE) target_node is the node that 1_best_candidate_node will be attached to. 

target_layer = target_node.Property8;

%% make one best candidate node
% check whether 1-best candidate node can be attached to the tree (Prperty5: number in layer)
if ~isempty(target_node.Next) % if the next of target_node is not empty 
    num_in_layer_for_one_best = target_node.Next.Property5 + 1;  
else   % if the next of target_node is empty 
    num_in_layer_for_one_best = 1; 
end

%% target list assignment
target_alpha_list = [];
remainder = rem(k,L);

if (k - remainder) == target_layer % from remaining layer to final layer
    target_alpha_list = sorted_alphaR_list;
    [ node, flag_ok_to_make ] = make_one_best_candidate_node_core_for_Rstep( target_node, target_alpha_list, num_in_layer_for_one_best, node, max_index, max_upper_bound_layer, k, L, remainder );
else
    if target_layer == 0
        target_alpha_list = sorted_alphaL_list_for_layerL;
    else
        target_alpha_list = sorted_alphaL_list;
    end
    [ node, flag_ok_to_make ] = make_one_best_candidate_node_core_for_Lstep( target_node, target_alpha_list, num_in_layer_for_one_best, node, max_index, max_upper_bound_layer, k, L, sorted_alphaL_list_for_layerL, sorted_alphaR_list );
end







end

