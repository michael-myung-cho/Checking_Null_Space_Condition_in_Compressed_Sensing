function [ output_args ] = display_updated_node_info_only( best_candidate_node, node, flag_ok_to_make )
%DISPLAY_UPDATED_NODES Summary of this function goes here
%   Detailed explanation goes here

     disp('***************************');
     disp('updated nodes are displayed');
     disp('***************************');
     dispDetail(best_candidate_node);
     disp('----------------------------');
     if flag_ok_to_make == 1 
        dispDetail(node(size(node,2)));             
     end
     
end

