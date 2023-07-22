classdef dlnode < hgsetget
    %LINK Summary of this class goes here
    %   Detailed explanation goes here
    properties
        Property0           % simple arithmetic
        Property1 = -1      % cheap bound
%        Property2 = -1      % linear_prog_algo 
        Property3 = -1      % alphaK
        Property4           % max(p1,p2,p3)
        Property5 = 1       % number in the layer 
        Property6 = 0       % node status '0':active, '1':submerged', '2':ignored
        Property7 = []      % index number
        Property8           % layer number
        Property9           % node number
        Next 
        Prev
    end
    methods
        function node = dlnode(Index)
            if nargin >= 0
                node.Property9 = Index;
            end
        end

        function add(newNode, beforeNode)
            newNode.Prev = beforeNode;
            beforeNode.Next = newNode;
        end
        function nextNode = disp(obj)
                fprintf('[%d]' , obj.Property7);
                nextNode = obj.Prev;
        end
        function dispDetail(obj)
            fprintf('p0: %d \n', obj.Property0);
            fprintf('p1: %d \n', obj.Property1);
%            fprintf('p2: %d \n', obj.Property2);
            fprintf('p3: %d \n', obj.Property3);
            fprintf('p4: %d \n', obj.Property4);
            fprintf('p5: %d \n', obj.Property5);
            fprintf('p6: %d \n', obj.Property6);
            fprintf('index: [%d] \n', obj.Property7);
            fprintf('layer level: %d \n', obj.Property8);            
            fprintf('num: %d \n', obj.Property9);
            obj.Prev
            obj.Next
            fprintf('\n');            
        end
    end
    
end


