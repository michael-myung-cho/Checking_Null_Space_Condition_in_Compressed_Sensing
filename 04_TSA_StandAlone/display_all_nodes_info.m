function [ output_args ] = display_all_nodes_info( node )
%DISPLAY_ALL_NODES Summary of this function goes here
%   Detailed explanation goes here

     disp('***************************');
     disp('All nodes are displayed');
     disp('***************************');
     for index_temp = 1: size(node,2)
         dispDetail(node(index_temp));
         disp('----------------------------');
     end

end

