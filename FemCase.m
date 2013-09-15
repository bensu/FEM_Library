classdef FemCase < hgsetget
    properties
        mesh
        SElist
        bc
        loads
        displacements
        reactions
        ncambio
        StressArray
        StrainArray
        StressArrayEle 
        StrainArrayEle
        MaxStress
        MaxStressInd
        MaxStressPoint
        MinStress
        MinStressInd
        MinStressPoint
        MaxStressEle
        MaxStressIndEle
        MaxStressPointEle
        MinStressEle
        MinStressIndEle
        MinStressPointEle
    end
    
    methods (Static)
   
        function [global_coordinates, ncambio] = nlocal2global(coordinates,coordinates2)
            ncambio = zeros(size(coordinates2,1),1);
            tol = 1e-10;
            maxnodes = size(coordinates,1);
            count = 1;
            global_coordinates = coordinates;
            for i = 1:size(coordinates2,1)
                for j = 1:maxnodes
                    if norm(coordinates2(i,:) - coordinates(j,:)) < tol
                        ncambio(i) = j;
                        break
                    end
                end
                if j == maxnodes
                    ncambio(i) = maxnodes + count;
                    global_coordinates = [global_coordinates; coordinates2(i,:)];
                    count = count + 1;
                end
            end
        end
        function global_connections = global_connect(connections,ncambio)
            global_connections = zeros(size(connections));
            for i = 1:size(connections,1)
                for j = 1:size(connections,2)
                    global_connections(i,j) = ncambio(connections(i,j));
                end
            end
        end
        function obj =  case_from_SElist(SElist)
            SE1 = SElist(1);
            SE2 = SElist(2);
            [global_coordinates, ncambio] = FemCase.nlocal2global(SE1.get('coordinates'),SE2.get('coordinates'));
            global_connections = FemCase.global_connect(SE1.get('connections'),ncambio);
            mesh_fem = MeshClass(global_coordinates,[SE1.get('connections');global_connections],SE1.get('material'));
            mesh_fem.set('ncambio',ncambio);
            obj = FemCase(mesh_fem,mesh_fem.sdof(),mesh_fem.sdof());
            SE1.set('ncambio',1:SE1.nnodes())
            SE2.set('ncambio',ncambio)
            obj.set('SElist',SElist)
        end
        function l_index1 = l_index(bsindex,index_out)
            l_index1 = zeros(length(index_out),1);
            for i = 1:length(index_out)
                l_index1(i) = find(bsindex == index_out(i));
            end
            
        end    
    end
    
    methods
        
        %% Initialize
        
        function obj = FemCase(meshin,bcin,loadsin)
            set(obj,'mesh',meshin);
            if strcmp(class(bcin),'BC')
                Bc = bcin;
            elseif isempty(bcin)
                Bc = BC(true(meshin.nnodes(),3));
            else Bc = BC(bcin);
            end
            set(obj,'bc',Bc);
            if isempty(loadsin)
                loadsin = zeros(meshin.nnodes(),3);
            else
            end
            set(obj,'loads',Loads(loadsin));
        end
        %% Solve

        function solve_condensed(obj)
            
            SElist = obj.get('SElist');
            SE1 = SElist(1);
            SE2 = SElist(2);
            
            index_in = SE1.nodesindex(SE1.nodes_in());
            index_out = SE1.nodesindex(SE1.nodes_out());
            index_in2 = SE2.nodesindex(SE2.nodes_in());
            index_out2 = SE2.nodesindex(SE2.nodes_out());
            
            stiffness = SE1.stiff();
            Kcc = stiffness(index_in,index_in);
            Krr = stiffness(index_out,index_out);
            Kcr = stiffness(index_in,index_out);
            Krc = Kcr';
            
            K_condensed = Krr- Krc*(Kcc\Kcr);
            
            nodos_unicos = unique([SE1.nodes_out() SE2.nodes_out()]);
            bsindex = SE1.nodesindex(nodos_unicos); %todos los grados de libertad del condensado.
            nodos_cond1 = FemCase.l_index(nodos_unicos,SE1.nodes_out());  %los nodos del 1er elemento respecto al condensado
            nodos_cond2 = FemCase.l_index(nodos_unicos,SE2.nodes_out());  % es decir, tiene indices mas chiquitos.
            
            l_index1 = SE1.nodesindex(nodos_cond1); %grados de libertd de esos indeces chicos.
            l_index2 = SE2.nodesindex(nodos_cond2);
            
            condensed_sdof = length(bsindex);
            
            SmallStiff = sparse(condensed_sdof,condensed_sdof);
            
            SmallStiff(l_index1,l_index1) = SmallStiff(l_index1,l_index1) + K_condensed;
            SmallStiff(l_index2,l_index2) = SmallStiff(l_index2,l_index2) + K_condensed;
            
            condensed_free = obj.get('bc').get('nodelist');
            condensed_free = condensed_free(bsindex);
            
            %             condensed_loadsv = reshape(loads',[],1);
            condensed_loadsv = obj.get('loads').get('nodelist');
            condensed_loadsv = condensed_loadsv(bsindex);
            
            condensed_dis = sparse(condensed_sdof,1);
            condensed_dis(condensed_free) = SmallStiff(condensed_free,condensed_free) \ condensed_loadsv(condensed_free);
            condensed_reactions = SmallStiff*condensed_dis;
            
            sdof = obj.get('mesh').sdof();
            displacements = zeros(sdof,1);
            displacements(index_out) = condensed_dis(l_index1); %desplazamientos de caras del 1er elemento
            displacements(index_in) = -(Kcc\Kcr)*displacements(index_out);%desplazamientos interiores del 1er elemento
            displacements(index_out2) = condensed_dis(l_index2); %desplazamientos de caras del 2do elemento
            displacements(index_in2) = -(Kcc\Kcr)*displacements(index_out2);%desplazamientos interiores del 2do elemento
            
            reactions = zeros(size(displacements));
            reactions(bsindex) = condensed_reactions;
            
            obj.set('displacements',VectorField(displacements));
            obj.set('reactions',VectorField(reactions));
        end
        
        function solve(obj)
            meshnow = obj.get('mesh');
            stiffness = meshnow.stiff();
            dis = zeros(meshnow.sdof(),1);
            free = obj.get('bc').get('nodelist');
            loadsv = obj.get('loads').get('nodelist');
            dis(free) = stiffness(free,free)\loadsv(free);
            dis = VectorField(dis);
            obj.set('displacements',dis);
            R = stiffness*dis.get('nodelist');
            obj.set('reactions',VectorField(R));
            %some sanity checks
            tol = 1e-15;
            for i = 1:3
                %sum(R(i:3:end))
                if sum(R(i:3:end))>tol
                    warning('sumR does not check')
                    sum(R(i:3:end))
                end
            end
        end
        %% Stress
        function tensors(obj)
            displacements = obj.get('displacements').get('nodelist');
            StressArray = zeros(obj.get('mesh').nnel(),8,6);
            StrainArray = StressArray;
            StressArrayEle = zeros(obj.get('mesh').nnel(),6);
            StrainArrayEle = StressArrayEle;
            gaussn = 2;
            [gaussp gaussw] = lgwt(gaussn,-1,1);
            for ele = 1:obj.get('mesh').nnel()
                new_element = obj.get('mesh').elementcreate(ele);
                count = 1;
                for i = 1:gaussn
                    xi = gaussp(i);
                    for j = 1:gaussn
                        eta = gaussp(j);
                        for k = 1:gaussn
                            mu = gaussp(k);
                            aux = new_element.Bmatrix(xi,eta,mu)*displacements(new_element.ldofsid());
                            C = new_element.get('material').get('C');
                            StrainArray(ele,count,:) = aux;
                            StressArray(ele,count,:) = C*aux;
                            count = count + 1;
                            StrainArrayEle(ele,:) = StrainArrayEle(ele,:) + aux'/8;
                            StressArrayEle(ele,:) = StressArrayEle(ele,:) + (C*aux)'/8;
                        end
                    end
                end
            end
            obj.set('StrainArray',StrainArray);
            obj.set('StressArray',StressArray);
            obj.set('StrainArrayEle',StrainArrayEle);
            obj.set('StressArrayEle',StressArrayEle);
        end
        function ExtremeStress(obj,opts)
            StressArrayEle = obj.get('StressArrayEle');
            StrainArrayEle = obj.get('StrainArrayEle');
            coordinates = obj.get('mesh').get('coordinates');
            MaxStress = zeros(3,1);
            MaxStressInd = zeros(3,1);
            MaxStressPoint = zeros(3,3);
            for i = 1:3
                switch opts
                    case 1
                        [val ind] = max(StressArray(:,:,i));
                        [val ind2] = max(val);
                    case 2
                        [val ind] = min(StressArray(:,:,i));
                        [val ind2] = min(val);
                end     
                MaxStress(i) = val;
                MaxStressInd(i) = ind(ind2);
                MaxStressPoint(i,:) = coordinates(connections(MaxStressInd(i),1),:); %Devuelve la coordenada del nodo uno de ese elemento
            end
            switch opts
                case 1
                    obj.set('MaxStress',MaxStress)
                    obj.set('MaxStressInd',MaxStressInd)
                    obj.set('MaxStressPoint',MaxStressPoint)
                case 2
                    obj.set('MinStress',MaxStressEle)
                    obj.set('MinStressInd',MaxStressInd)
                    obj.set('MinStressPoint',MaxStressPoint)
            end
        end
        function ExtremeStressElement(obj,opts)
            StressArrayEle = obj.get('StressArrayEle');
            StrainArrayEle = obj.get('StrainArrayEle');
            coordinates = obj.get('mesh').get('coordinates');
            connections = obj.get('mesh').get('connections');
            MaxStressEle = zeros(3,1);
            MaxStressIndEle = zeros(3,1);
            MaxStressPointEle = zeros(3,3);
            for i = 1:3
                switch opts
                    case 1
                        [val ind] = max(StressArrayEle(:,i));
                    case 2
                        [val ind] = min(StressArrayEle(:,i));
                end   
                MaxStressEle(i) = val;
                MaxStressIndEle(i) = ind;
                MaxStressPointEle(i,:) = coordinates(connections(ind,1),:); %Devuelve la coordenada del nodo uno de ese elemento
            end
            switch opts
                case 1
                    obj.set('MaxStressEle',MaxStressEle)
                    obj.set('MaxStressIndEle',MaxStressIndEle)
                    obj.set('MaxStressPointEle',MaxStressPointEle)
                case 2
                    obj.set('MinStressEle',MaxStressEle)
                    obj.set('MinStressIndEle',MaxStressIndEle)
                    obj.set('MinStressPointEle',MaxStressPointEle)
            end
        end
        
        function getMaxStressEle(obj)
            obj.tensors();
            obj.ExtremeStressElement(1)
        end
            
        %% Plot
        
        function plot(obj,elelist)
            scale = 1000;
            hold on
            mesh = obj.get('mesh');
            mesh = MeshClass(mesh.get('coordinates'),mesh.get('connections'),mesh.get('material'));
            dis = obj.get('displacements');
            Dis = dis.xyzout();
            mesh.plot('g');
            bc = obj.get('bc');
            bc = VectorField(~bc.get('nodelist'));
            bc.plotVF(mesh.get('coordinates'),'k')
            if ~isempty(dis)
                dis.plotVF(mesh.get('coordinates'),'b')
            end
            obj.get('loads').plotVF(mesh.get('coordinates'),'r');
            hold off
        end
        
        
    end
    
    
end
