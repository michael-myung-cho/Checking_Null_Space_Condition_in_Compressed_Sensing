function [ sorted_table_for_p4_tracking ] = sorting_operation_from_bottom( table_for_p4_tracking )
%SORTING_OPERATION Summary of this function goes here
%   Detailed explanation goes here

temp = table_for_p4_tracking(end,:);
table_for_p4_tracking(end,:)=[];
[a,b] = size(table_for_p4_tracking);
if a == 0 % table size = 1
    sorted_table_for_p4_tracking = temp;
elseif a == 1
    if table_for_p4_tracking(1,2) < temp(1,2)
        sorted_table_for_p4_tracking = [table_for_p4_tracking;temp];
    else
        sorted_table_for_p4_tracking = [temp;table_for_p4_tracking];
    end
elseif a >= 2
    % ex)       (updated) 5| 3| 4| 5| 7| 8| 11| 15 ordered
    % actual table index: 1| 2| 3| 4| 5| 6|  7|  8
    % used table index  : x| 1| 2| 3| 4| 5|  6|  7
    % actual table index 1 is treated as target index. Thus, it is not
    % used for table index
    target_value = temp(1,2);
    first_index = 1;
    last_index = a;
    while (first_index < last_index)
        mid_index = floor((first_index + last_index)/2);
        mid_value = table_for_p4_tracking(mid_index,2);
        if (target_value < mid_value)
            last_index = mid_index;
        else
            first_index = mid_index + 1;
        end 
    end
    if first_index == 1
            sorted_table_for_p4_tracking = [temp;table_for_p4_tracking];       
    elseif first_index == a
            if table_for_p4_tracking(a,2) < temp(1,2)
                sorted_table_for_p4_tracking = [table_for_p4_tracking;temp];
            else
                sorted_table_for_p4_tracking = [table_for_p4_tracking(1:end-1,:);temp;table_for_p4_tracking(end,:)];
            end
    else
        if table_for_p4_tracking(first_index-1,2) <= temp(1,2) && table_for_p4_tracking(first_index,2) >= temp(1,2)
            sorted_table_for_p4_tracking = [table_for_p4_tracking(1:first_index-1,:);temp; table_for_p4_tracking(first_index:end,:)];
        elseif table_for_p4_tracking(first_index,2) <= temp(1,2) && table_for_p4_tracking(first_index+1,2) >= temp(1,2)
            sorted_table_for_p4_tracking = [table_for_p4_tracking(1:first_index,:);temp; table_for_p4_tracking(first_index+1:end,:)];
        else
        end
    end
end
end

