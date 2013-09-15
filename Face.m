classdef Face < hgsetget
    %A face is described by a list of element faces
    %each element face is described by it's element id, and the plane
    %coordnum == value
    properties 
        nodelist
        parent
        elementlist
        coordnumlist
        valuelist
        par_plot
    end
    methods (Static)
        function newface = new_face(parent,coordnum,value)
            tol = 1e-15;
            coords = parent.get('coordinates');
            nodelist = [];
            for nn = 1:parent.nnodes()
                if abs(coords(nn,coordnum)-value)<tol
                    nodelist = [nodelist;nn];
                end
            end
            if or(coordnum > 3,coordnum<1)
                error(strcat('Invalid coordnumber:',num2str(coordnum)))
            end
            newface = Face(nodelist,parent);
        end
       
    end
    
    methods
        function obj = Face(nodelist,parent)
            if ~strcmp('MeshClass',class(parent))
                error('wrong type for parent')
            end
            set(obj,'nodelist',nodelist);
            set(obj,'parent',parent);
            parent.add_outer_face(obj)
        end
        
        function [facevalue facecoord] = normalnodes2valcoord(obj,normalnodes)
                M = max(normalnodes);
                if M==8
                    aux = 1;
                else aux = -1;
                end
                facevalue = aux;
                switch M-min(normalnodes)
                    case 6
                        aux = 1;
                    case 5
                        aux = 2;
                    case 3
                        aux = 3;
                end
                facecoord = aux;
        end
        
        function [elelist facecoord facevalue] = faceinelement(obj)
            parent1 = obj.get('parent');
            connections = parent1.get('connections');
            nface = obj.get('nodelist');
            elematrix = zeros(size(connections,1),4); %MAGIC NUMBER 4
            count = ones(size(connections,1),1);
            for i = 1:length(nface)
                elenodes = parent1.elementsofnode(nface(i));
                for j = 1:length(elenodes)
                    el = elenodes(j);
                    elematrix(el,count(el)) = nface(i);
                end
                count(elenodes) = count(elenodes) + 1;
            end
            
            %loop to transform get out of the elements the value and coordinates        
            elelist = [];  
            for i = 1:size(elematrix,1)
                if all(elematrix(i,:))
                    elelist = [elelist i];
                else elematrix(i,:) = 0;
                end
            end
            facevalue = zeros(size(elelist));
            facecoord = facevalue;
            
            for i = 1:length(elelist);
                elenum = elelist(i);
                foundnodes = elematrix(elenum,:);   %should be 4x1
                allnodes = connections(elenum,:) ;  %should be 8x1
                normalnodes = zeros(size(foundnodes));
                for j = 1:length(foundnodes)
                    normalnodes(j) = find(foundnodes(j) == allnodes);
                end

                [facevalue_aux facecoord_aux] = normalnodes2valcoord(obj,normalnodes);
                facevalue(i) = facevalue_aux;
                facecoord(i) = facecoord_aux;

            end
            
            obj.set('elementlist',elelist);
            obj.set('coordnumlist',facecoord);
            obj.set('valuelist',facevalue);
        end  
    end
end