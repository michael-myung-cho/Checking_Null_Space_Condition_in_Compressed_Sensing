function [ table_for_p4_tracking ] = update_table_from_node_new( table_for_p4_tracking, best_candidate_node, flag_ok_to_make, node )
% UPDATE_TABLE_FROM_THE_NODE_NEW
% tracking the best candidate node

if ((best_candidate_node.Property6 == 1) || (best_candidate_node.Property6 == 2))
    table_for_p4_tracking(1,:) = []; % remove from table
else
    table_for_p4_tracking(1,2) = best_candidate_node.Property4; % update
    % sorting operation from top means conducting sorting operation for
    % best candidate node which has been updated in p4
    table_for_p4_tracking = sorting_operation_from_top( table_for_p4_tracking ); % top means 1st row of the table
end
if flag_ok_to_make == 0 % no node is attached in the tree
     % nothing happen
else % one node is attached in the tree
    [a,b] = size(table_for_p4_tracking);
    one_best_candidate_node = node(size(node,2));
    table_for_p4_tracking(a+1,1) = one_best_candidate_node.Property9;
    table_for_p4_tracking(a+1,2) = one_best_candidate_node.Property4;
    % sorting operation from bottom means conducting sorting operation for
    % a node newly attached to the tree
    table_for_p4_tracking = sorting_operation_from_bottom( table_for_p4_tracking ); % bottom means the last row of the table
end
end

