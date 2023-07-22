function [ output_args ] = display_tree( node )
%DISPLAY_TREE Summary of this function goes here
%   Detailed explanation goes here

nums_of_nodes = size(node,2);

% find max. layer (level)
level = 0;
level_nodes_temp = zeros(nums_of_nodes,1);

for i_disp_nodes = 1:nums_of_nodes
    level_nodes_temp(i_disp_nodes,1) = node(i_disp_nodes).Property8; % layer
end
level = max(level_nodes_temp);

disp('tree display ===============================');

while (level ~= 0)
    
    for i = 1:nums_of_nodes
        if (isempty(node(i).Next)) && (node(i).Property8 == level)

            current_node = node(i);
        
            for p = level:-1:1
                
                if current_node.Property8 == p
                    next_node = disp(current_node);
                    fprintf(' -> ');
                    current_node = next_node;
                else
                    fprintf(' -> ');                   
                end       

            end
            fprintf(' root '); % root node display
            fprintf('\n');
        end
    end
    
    level = level - 1;
        
end
disp('============================================');
    
end


