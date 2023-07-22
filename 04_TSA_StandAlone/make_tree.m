function [ alphaTSA_LB, alphaTSA_UB ] = make_tree( permutated_A, sorted_alphaL_list, sorted_alphaR_list, sorted_alphaL_list_for_layerL, kv, L)

[m,n] = size(permutated_A);
%% save data to xls file
% iCT = 1; % current count
% clk=clock;
% datestmp=strcat([date,'-',num2str(clk(4)),'h',num2str(clk(5))]);
% filenameTSA = strcat('(',num2str(m),'x',num2str(n),')','TSA_Mid_Result-',datestmp,'.xlsx');
% col_name = {'m,','n','kv','itr','alphaTSA_LB','alphaTSA_UB';};
% xlswrite(filenameTSA,col_name);
%% time
maxTime=60*60*24; % 1day

%% determine sparsity
max_index = size(permutated_A,2);
permutated_H = null(permutated_A);
%% root node assignment (initial best candidate node)
initial_vector = zeros(1,L);
node(1) = dlnode(1);   % root number = 1
node(1).Property1 = 0; % cheap upper bound
%node(1).Property2 = 0; % linear upper bound
node(1).Property3 = 0; % exact alphaK
node(1).Property4 = 0; % main metric
node(1).Property6 = 0; % status: submerged
node(1).Property7 = initial_vector; % root node index = 0
node(1).Property8 = 0; % layer 0

%% table for lowerest property4 tracking
table_for_p4_tracking = zeros(1,2);
table_for_p4_tracking(1,1) = 1; % node number
table_for_p4_tracking(1,2) = node(1).Property4; 
currnet_smallest_p4_node_number = 0;

%% calculate maximum upper bound on layer-k
% max_upper_bound_layer = sum(sorted_alphaL_list_for_layerL(1:floor(k/L),end)) + sum(sorted_alphaR_list(1:rem(k,L),end));
max_upper_bound_layer = 0;
alphaTSA_LB = 0;
flagF = 0;
%% Step 3
tTSA = tic;
itr=1;
while (1) 
    flag_ok_to_make = 0;
    
    %% find best candidate node
    best_candidate_node = node(table_for_p4_tracking(1,1));
    num_layer = best_candidate_node.Property8;
%     disp('***************************');
%     disp('(Step2) best_candidate_node');
%     disp('***************************');
%    dispDetail(best_candidate_node);
    
    %% Step 3-1: making node and linking it to tree        
    if ( best_candidate_node.Property1 == -1) % if P1 NOT calculated (cheap bound operation)
        step = 31;
%         disp('***************************');
%         disp('Step 3-1');
%         disp('***************************');
                
        % cheap bound calculation instead of linear programming due to operation time reduction
        best_candidate_node.Property1 = cheap_bound_algo(best_candidate_node, sorted_alphaL_list_for_layerL, L);        
        
        % calculate remaining upper bound over k-j elements
        sum_remaining_alpha = remaining_upper_bound_BCN( num_layer, kv, L, best_candidate_node, sorted_alphaL_list_for_layerL, sorted_alphaR_list );
        
        % find minimum property value of the parent node
        min_best_node = best_candidate_node.Property0;
        if (best_candidate_node.Property1 ~= -1) && (min_best_node > best_candidate_node.Property1)
             min_best_node = best_candidate_node.Property1;
        end
%         if (best_candidate_node.Property2 ~= -1) && (min_best_node > best_candidate_node.Property2)
%              min_best_node = best_candidate_node.Property2;
%         end
        if (best_candidate_node.Property3 ~= -1) && (min_best_node > best_candidate_node.Property3)
             min_best_node = best_candidate_node.Property3;
        end
        
        % update P4 of best candidate node
        best_candidate_node.Property4 = max_upper_bound_layer - (min_best_node + sum_remaining_alpha);

        % make 1-best-candidate-node in parent node
        % (NOTE) target_node is the node that 1_best_candidate_node will be attached to.
        if (num_layer == 0)
        else
%            fprintf('node is attached to the parent node\n');
            target_node = best_candidate_node.Prev;
            [node, flag_ok_to_make] = make_one_best_candidate_node( target_node, node, sorted_alphaL_list, sorted_alphaR_list, sorted_alphaL_list_for_layerL, kv, L, max_index, max_upper_bound_layer );
        end
  
%     %% Step 3-2: calculate p2 (linear programming)  
%     elseif (best_candidate_node.Property2 == -1) %if P2 calcuated but NOT P3
%         step = 32;
%         disp('***************************');
%         disp('Step 3-2');
%         disp('***************************');
%         if (num_layer == kv)   
%             % linear programming algorithm for obtaining tighter bound
% %            linear_alphaK =  linear_prog_algo(sorted_alpha1_list, sorted_alpha2_list, best_candidate_node);
% %            best_candidate_node.Property2 = linear_alphaK; % calculate p2
%             best_candidate_node.Property2 = best_candidate_node.Property1; %temp
% 
%             % calculate remaining upper bound over k-j elements
%             sum_remaining_alpha = remaining_upper_bound( num_layer, kv, L, best_candidate_node, sorted_alphaL_list_for_layerL, sorted_alphaR_list );
%             
%             % find minimum property value of the parent node
%             min_best_node = best_candidate_node.Property0;
%             if (best_candidate_node.Property1 ~= -1) && (min_best_node > best_candidate_node.Property1)
%                  min_best_node = best_candidate_node.Property1;
%             end
%             if (best_candidate_node.Property2 ~= -1) && (min_best_node > best_candidate_node.Property2)
%                  min_best_node = best_candidate_node.Property2;
%             end
%             if (best_candidate_node.Property3 ~= -1) && (min_best_node > best_candidate_node.Property3)
%                 min_best_node = best_candidate_node.Property3;
%             end
%             
%             % update P4 of best candidate node
%             best_candidate_node.Property4 = max_upper_bound_layer - (min_best_node + sum_remaining_alpha);       
%         else
%             % not on final layer
%             best_candidate_node.Property2 = best_candidate_node.Property1;
%         end
        
    %% Step 3-3: make 1_best_candidate_node    
    elseif (best_candidate_node.Property3 == -1 ) %if P2 calcuated but NOT P3
        step = 33;
%         disp('***************************');
%         disp('Step 3-3');
%         disp('***************************');
        best_candidate_node.Property3 = exact_alphaK_algo( permutated_A, permutated_H, sorted_alphaL_list, best_candidate_node, kv, L );
        % exact alphaK is the smallest value among all upper bound.
        min_best_node = best_candidate_node.Property3;
        % calculate remaining upper bound over k-j elements
        sum_remaining_alpha = remaining_upper_bound_BCN( num_layer, kv, L, best_candidate_node, sorted_alphaL_list_for_layerL, sorted_alphaR_list );          
        % update P4 of best candidate node: it is on the final layer,
        % hence, it does not have remaining part
        best_candidate_node.Property4 = max_upper_bound_layer - (min_best_node + sum_remaining_alpha); 
        
        if (num_layer == kv) && alphaTSA_LB < min_best_node
            alphaTSA_LB = min_best_node;
        end
    %% Step 3-4: When P3 has already been calculated ("submerge a node" operation)
    else
        best_candidate_node.Property6 = 1; % status: submerged
        if (num_layer == kv)
            currnet_smallest_p4_node_number = best_candidate_node.Property9; % num of node
            alphaTSA_UB = node(currnet_smallest_p4_node_number).Property3;
            break;
        else
            % make 1-best-candidate-node in current node(best-candidate-node)
            % (NOTE) target_node is the node that 1_best_candidate_node will be attached to.
%            fprintf('node is attached to the current node\n');            
            target_node = best_candidate_node;
            [node, flag_ok_to_make] =  make_one_best_candidate_node( target_node, node, sorted_alphaL_list, sorted_alphaR_list, sorted_alphaL_list_for_layerL, kv, L, max_index, max_upper_bound_layer );  
        end
        

    end
  
    %% display updated nodes only
%    display_all_nodes_info(node);
%    display_updated_node_info_only(best_candidate_node, node, flag_ok_to_make);

    %% tree display
%    display_tree(node);    
    

    %% update table for smallest p4     
    % new sorting operation is implemented
    table_for_p4_tracking = update_table_from_node_new(table_for_p4_tracking, best_candidate_node, flag_ok_to_make, node);

    t1TSA = toc(tTSA);    
    alphaTSA_UB = -table_for_p4_tracking(1,2);
    if(mod(itr,500)==1)
        fprintf('TSA time(min.): %d (LB = %f, UB = %f)\n',t1TSA,alphaTSA_LB,alphaTSA_UB);
    end
    itr=itr+1;
    %% time break
    if t1TSA > maxTime
        break;
    end
    %% intermidiate result save
%    if (t1TSA >= unitT*iCT)
%     if(mod(itr,500)==1)
%         iCT = iCT + 1;
%         % save data to xls file
%         resbufTSA=[m, n, kv, itr, alphaTSA_LB, alphaTSA_UB];
%         xlRangeTSAM = sprintf('A%d',iCT+1);
%         xlswrite(filenameTSA,resbufTSA,'Sheet1',xlRangeTSAM);
%     end
%    itr=itr+1;
%    end
%     if (t1TSA >= tT)
%         break;
%     end
end
fprintf('TSA time(min.): %d (LB = %f, UB = %f)\n',t1TSA,alphaTSA_LB,alphaTSA_UB);


%% save data to xls file
% iCT = iCT + 1;
% resbufTSA=[m, n, kv, alphaTSA_LB, alphaTSA_UB, t1TSA];
% xlRangeTSAM = sprintf('A%d',iCT+1);
% xlswrite(filenameTSA,resbufTSA,'Sheet1',xlRangeTSAM);


%% display node property
% display_all_nodes_info(node);
% display_tree(node); 

%% display information
%  disp('K set node is found.');
%  dispDetail( node(currnet_smallest_p4_node_number) );
%  current_node = node(currnet_smallest_p4_node_number);    
%  for i_node = 1:node(currnet_smallest_p4_node_number).Property8
%      if  current_node.Property8 ~= 0 
%          next_node = disp(current_node);
%          fprintf(' -> ' );
%          current_node = next_node;
%      end
%  end
%  fprintf('\n');
%  fprintf('Total number of nodes: %d\n', size(node,2));
%  diary off
%  save('TSA_3step_all_info.mat');
end


