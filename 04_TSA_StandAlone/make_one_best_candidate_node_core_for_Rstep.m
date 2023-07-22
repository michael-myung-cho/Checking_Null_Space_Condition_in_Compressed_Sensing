function [ node, flag_ok_to_make ] = make_one_best_candidate_node_core_for_Rstep( target_node, target_alpha_list, num_in_layer_for_one_best, node, max_index, max_upper_bound_layer, k, L, remainder )
% MAKE_ONE_BEST_CANDIDATE_NODE_CORE_FOR_RSTEP \
% making one_best_candidate_node on remaining layer to reach the final layer

target_layer = target_node.Property8;
flag_ok_to_make = 0;

while(1)
    if num_in_layer_for_one_best <= size(target_alpha_list,1)
        % check whether tree order is not violated or not
        if ( abs(target_node.Property7(1,end)) + k - target_layer <= max_index ...
             && abs(target_alpha_list(num_in_layer_for_one_best,1)) > abs(target_node.Property7(1,end)) ) % max(num_in_layer_for_one_best) = max_index

            flag_ok_to_make = 1;            
            num_layer_for_one_best = target_layer + remainder; % because of remaining step
 
            nums_of_nodes = size(node,2);
            node(nums_of_nodes+1) = dlnode(nums_of_nodes+1);
            one_best_candidate_node = node(nums_of_nodes+1);
            one_best_candidate_node.Prev = target_node;
            target_node.Next = one_best_candidate_node;
            
            % in order to make index [2,3,3] in 3 step tree search with
            % 2 remaining step
            index_for_remaining = [];
            for i_in = 1:L-remainder
                index_for_remaining = [target_alpha_list(num_in_layer_for_one_best, end-1), index_for_remaining];
            end
            
            one_best_candidate_node.Property7 = [target_alpha_list(num_in_layer_for_one_best,1:end-1), index_for_remaining];
            one_best_candidate_node.Property8 = num_layer_for_one_best;
            one_best_candidate_node.Property5 = num_in_layer_for_one_best;
                          

            % find minimum property value of the parent node
            min_parent_node = target_node.Property0;
            if (target_node.Property1 ~= -1) && (min_parent_node > target_node.Property1)
                min_parent_node = target_node.Property1;
            end
%             if (target_node.Property2 ~= -1) && (min_parent_node > target_node.Property2)
%                 min_parent_node = target_node.Property2;
%             end
            if (target_node.Property3 ~= -1) && (min_parent_node > target_node.Property3)
                min_parent_node = target_node.Property3;
            end
                
            % P0 = min(P0,P1,P2,P3) of parent node + alpha_{l,L}  
            % one best candidate node is on final layer, hence, it does not
            % have remaining part to be considered.
            one_best_candidate_node.Property0 = min_parent_node + target_alpha_list(num_in_layer_for_one_best, end);
            one_best_candidate_node.Property4 = max_upper_bound_layer -(one_best_candidate_node.Property0);
            break;
        else
            num_in_layer_for_one_best = num_in_layer_for_one_best + 1;
        end
    else
       target_node.Property6 = 2; % status ignored
       break;
    end
end

end

