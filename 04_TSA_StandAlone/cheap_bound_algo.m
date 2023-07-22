function [ output_cheap_bound ] = cheap_bound_algo( current_node, sorted_alphaL_list, L )
% cheap_bound_algo - calculate cheap upper bound for given indices

current_layer = current_node.Property8; % num_layer = number of index
index_for_cheap_bound = [];

%% calculate cheap bound from given indexs
if current_layer == 0 % root node
    output_cheap_bound = 0;
else
    % find indices
    number_of_layers = ceil(current_layer/L);
    for i_in = 1:number_of_layers
        index_for_cheap_bound = [current_node.Property7, index_for_cheap_bound];
        current_node = current_node.Prev;
    end
    
    actual_index = index_for_cheap_bound(1,1:current_layer);        
    IN_L = nchoosek( actual_index, L );
            
    k_set_sum = 0;
    for q = 1:size(IN_L,1)  % IN_L(x,1) < IN2(x,2)
        if sign(IN_L(q,1)) == 1
            pickL_portion = find_portion(sorted_alphaL_list, IN_L(q,:),L);
        else
            pickL_portion = find_portion(sorted_alphaL_list, -IN_L(q,:),L);
        end
        k_set_sum = k_set_sum + pickL_portion;
    end
    output_cheap_bound = k_set_sum / nchoosek(size(actual_index,2)-1, L-1);  
end

end

