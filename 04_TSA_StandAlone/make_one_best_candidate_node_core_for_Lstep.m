function [ node, flag_ok_to_make ] = make_one_best_candidate_node_core_for_Lstep( target_node, target_alpha_list, num_in_layer_for_one_best, node, max_index, max_upper_bound_layer, k, L, sorted_alphaL_list_for_layerL, sorted_alphaR_list)
% MAKE_ONE_BEST_CANDIDATE_NODE_CORE_FOR_LSTEP 
% making one_best_candidate_node on regular layer
% (NOTE) target_node is the node that 1_best_candidate_node will be attached to.

    target_layer = target_node.Property8;
    flag_ok_to_make = 0;
    while(1)
        if num_in_layer_for_one_best < size(target_alpha_list,1)
            % check whether tree order is not violated or not
            if (( abs(target_node.Property7(1,end)) < abs(target_alpha_list(num_in_layer_for_one_best,1))) ...
                && (( abs(target_alpha_list(num_in_layer_for_one_best,end-1)) + k - (target_layer+L)) <= max_index) )

                flag_ok_to_make = 1;            
                num_layer_for_one_best = target_layer + L; % because of L step
 
                nums_of_nodes = size(node,2);
                node(nums_of_nodes+1) = dlnode(nums_of_nodes+1);
                one_best_candidate_node = node(nums_of_nodes+1);
                one_best_candidate_node.Prev = target_node;
                target_node.Next = one_best_candidate_node;
                one_best_candidate_node.Property7 = target_alpha_list(num_in_layer_for_one_best,1:end-1);
                one_best_candidate_node.Property8 = num_layer_for_one_best;
                one_best_candidate_node.Property5 = num_in_layer_for_one_best;              

                % calculate remaining upper bound over k-j elements
                sum_remaining_alpha = remaining_upper_bound_1BCN( num_layer_for_one_best, k, L, one_best_candidate_node, sorted_alphaL_list_for_layerL, sorted_alphaR_list );
                        
                % when one best candidate node is atteched in layerL (exact alphaL is already calculated) 
                if num_layer_for_one_best == L
                    exact_alphaL = target_alpha_list(num_in_layer_for_one_best,end);
                    one_best_candidate_node.Property0 = exact_alphaL;
                else
                    % find minimum property value of the parent node
                    min_parent_node = target_node.Property0;
                    if (target_node.Property1 ~= -1) && (min_parent_node > target_node.Property1)
                        min_parent_node = target_node.Property1;
                    end
%                     if (target_node.Property2 ~= -1) && (min_parent_node > target_node.Property2)
%                         min_parent_node = target_node.Property2;
%                     end
                    if (target_node.Property3 ~= -1) && (min_parent_node > target_node.Property3)
                        min_parent_node = target_node.Property3;
                    end
                    
                    % P0 = min(P0,P1,P2,P3) of parent node + alpha_{l,L}
                    exact_alphaL = target_alpha_list(num_in_layer_for_one_best,end) ;
                    one_best_candidate_node.Property0 = min_parent_node + exact_alphaL;
                end
                
                one_best_candidate_node.Property4 = max_upper_bound_layer -(one_best_candidate_node.Property0 + sum_remaining_alpha);
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

