%% Analisis Previo
clc
clear

save_on = true;
if save_on
% All units are in N, m.

tometer = 1e3;
b = 50/tometer;
t = 5/tometer;
H = 100/tometer;
L = 2*100/tometer;

h = (H-2*t)/2;

Inertia = b*(t^3)/12 + b*t*(h^2) + t*((H-2*t)^3)/12;
Area = ((H-2*t)*t+2*b*t);

P = -100*1000; %en Newtons
E = 200e09;	% Youngs modulus in Pascal.
nu = 0.3;
Dprevio = P*(L^3)/(3*E*Inertia); %Maximum Displacement, z-coodinate
phiprevio = P*(L^2)/(2*E*Inertia);
sigmaxx = P*L*(H/2)/Inertia;

%% Material

steel = Mat(1,'Iso',E,nu,1);

%% SE1

SE1 = MeshClass.mesh_import('coords.malladoNX.txt','connect.malladoNX.txt',steel);

coordnum = 2;   %in y axis.
bcface = Face.new_face(SE1,coordnum,SE1.mincoordval(coordnum));

coordnum = 2;   %in y axis.
middleface = Face.new_face(SE1,coordnum,SE1.maxcoordval(coordnum));

SE1.set('nodesout',[bcface.get('nodelist');middleface.get('nodelist')]);

%% SE2

SE2 = MeshClass.mesh_import('coords.malladoNX.txt','connect.malladoNX.txt',steel);

coordnum = 2;   %in y axis.
bcface = Face.new_face(SE2,coordnum,SE2.mincoordval(coordnum));
coordnum = 2;   %in y axis.
middleface = Face.new_face(SE1,coordnum,SE1.maxcoordval(coordnum));

SE2.set('nodesout',[bcface.get('nodelist');middleface.get('nodelist')]);
coordnum = 2;
SE2.shift([0 SE1.maxcoordval(coordnum) 0]);
case_NX = FemCase.case_from_SElist([SE1 SE2]);

%% BC
coordnum = 2;
bc_face = Face.new_face(case_NX.get('mesh'),coordnum,case_NX.get('mesh').mincoordval(coordnum));
case_NX.get('bc').addfixedface(bc_face);

%% Loads
coordnum = 2;
load_face = Face.new_face(case_NX.get('mesh'),coordnum,case_NX.get('mesh').maxcoordval(coordnum));
case_NX.get('loads').qinface(case_NX.get('mesh'),load_face,[0 0 P/Area]);

%% Solve
case_NX.solve_condensed()
case_NX.getMaxStressEle()

%% Save
save('case_NX')
else load 'case_NX.mat'
end

%% Plot

% 
% 
% % Stress
% 
% 
% 
% Maxx = case_NX.get('MaxStressEle');
% 
% SA = case_NX.get('StressArrayEle');
% case_NX.get('MaxStressIndEle')
% Dis = case_NX.get('displacements').xyzout();
% tita_out = (180/pi)*max(Dis(:,2))/(H/2);
% min(Dis(:,3))*tometer;
% 
% case_NX.get('MaxStressIndEle')
% 
% coordinates = case_NX.get('mesh').coordinates;
% 
% connections = case_NX.get('mesh').connections;
% 
% nod = case_NX.get('mesh').findnode([.0225 0 0.011]);
% 
% 
% case_NX.get('mesh').coordinates(case_NX.get('MaxStressIndEle'));
% 
% mesh = case_NX.get('mesh');
% 
% elelist = [];
% xx = unique(coordinates(:,3));
% for i = 1:length(xx)
%     x0 = [0.0225,0,xx(i)];
%     nodes = mesh.findnode(x0);
%     for j = 1:length(nodes)
%         elelist = [elelist mesh.elementsofnode(nodes)];
%     end
%     
%            
% end
% 
% elelist = unique(elelist);
% elecoords  = [];
% for i = 1:length(elelist)
%     elecoords = [elecoords; mesh.elementcoord(elelist(i))];
% end
% 
% 
% 
%             dis = case_NX.get('displacements');
%             Dis = dis.xyzout();
%             scale = 10;
%             hold on
%             mesh = case_NX.get('mesh');
%             mesh = MeshClass(mesh.get('coordinates')+Dis*scale,mesh.get('connections'),mesh.get('material'));
% 
%             case_NX.set('mesh',mesh')
%             case_NX.plot(case_NX.get('MaxStressIndEle'))
% plot(sort(elecoords(:,3)),sort(SA(elelist,2)))
% %
% elelist2 = [];
% xx = unique(coordinates(:,3));
% for i = 1:length(xx)
%     x0 = [0.0225,200/1e3,xx(i)];
%     nodes = mesh.findnode(x0);
%     for j = 1:length(nodes)
%         elelist2 = [elelist2 mesh.elementsofnode(nodes)];
%     end
%     
%            
% end
% 
% elelist2 = unique(elelist2);
% elecoords2  = [];
% for i = 1:length(elelist2)
%     elecoords2 = [elecoords2; mesh.elementcoord(elelist2(i))];
% end
% 
% %
% elelist3 = [];
% xx = unique(coordinates(:,3));
% for i = 1:length(xx)
%     x0 = [0.0225,100/1e3,xx(i)];
%     nodes = mesh.findnode(x0);
%     for j = 1:length(nodes)
%         elelist3 = [elelist3 mesh.elementsofnode(nodes)];
%     end
%     
%            
% end
% 
% elelist3 = unique(elelist3);
% elecoords3  = [];
% for i = 1:length(elelist3   )
%     elecoords3 = [elecoords3; mesh.elementcoord(elelist3(i))];
% end
% 
% %
% elelist4 = [];
% xx = unique(coordinates(:,3));
% for i = 1:length(xx)
%     x0 = [0.0225,50/1e3,xx(i)];
%     nodes = mesh.findnode(x0);
%     for j = 1:length(nodes)
%         elelist4 = [elelist4 mesh.elementsofnode(nodes)];
%     end
%     
%            
% end
% 
% elelist4 = unique(elelist4);
% elecoords4  = [];
% for i = 1:length(elelist4)
%     elecoords4 = [elecoords4; mesh.elementcoord(elelist4(i))];
% end
% 
% %
% elelist5 = [];
% xx = unique(coordinates(:,3));
% for i = 1:length(xx)
%     x0 = [0.0225,150/1e3,xx(i)];
%     nodes = mesh.findnode(x0);
%     for j = 1:length(nodes)
%         elelist5 = [elelist5 mesh.elementsofnode(nodes)];
%     end
%     
%            
% end
% 
% elelist5 = unique(elelist5);
% elecoords5  = [];
% for i = 1:length(elelist5)
%     elecoords5 = [elecoords5; mesh.elementcoord(elelist5(i))];
% end
% 
% plot(sort(elecoords(:,3)),sort(SA(elelist,2)),'o')
% 
% %
% 
% % case_NX.plot(case_NX.get('MaxStressIndEle'))
%             dis = case_NX.get('displacements');
%             Dis = dis.xyzout();
%             scale = 10;
%             hold on
%             mesh = case_NX.get('mesh');
%             mesh = MeshClass(mesh.get('coordinates')+Dis*scale,mesh.get('connections'),mesh.get('material'));
% 
%             mesh.plot('g');
%           elelist = case_NX.get('MaxStressIndEle');
%           elelist = elelist([2 3 1]);
%         cl = ['r','b','k'];
%           for i = 1:length(elelist)
%                 hold on
%                 element = elelist(i);
%                 coords = mesh.elementcreate(element).coordinates();
%                 elementcoord = mesh.elementcoord(element);
%                 quiver3(elementcoord(1),0*elementcoord(2),elementcoord(3),0,SA(element,2)/1e10,0,'r')
%           end
%           
%           for i = 1:length(elelist2)
%                 hold on
%                 element = elelist2(i);
%                 coords = mesh.elementcreate(element).coordinates();
%                 elementcoord = mesh.elementcoord(element);
%                 quiver3(elementcoord(1),200/1e3,elementcoord(3),0,SA(element,2)/1e10,0,'r')
%           end
%           max(SA(elelist,2))
%           for i = 1:length(elelist3)
%                 hold on
%                 element = elelist3(i);
%                 coords = mesh.elementcreate(element).coordinates();
%                 elementcoord = mesh.elementcoord(element);
%                 quiver3(elementcoord(1),100/1e3,elementcoord(3),0,SA(element,2)/1e10,0,'r')
%           end
%      
%           %
%           
%           case_NX.plot(case_NX.get('MaxStressIndEle'))
%             dis = case_NX.get('displacements');
%             Dis = dis.xyzout();
%             scale = 1untitled0;
%             hold on
%             mesh = case_NX.get('mesh');
%             mesh = MeshClass(mesh.get('coordinates')+Dis*scale,mesh.get('connections'),mesh.get('material'));
% 
%             mesh.plot('g');
%           for i = 1:length(elelist)
%                 hold on
%                 element = elelist(i);
%                 coords = mesh.elementcreate(element).coordinates();
%                 elementcoord = mesh.elementcoord(element);
%                 quiver3(elementcoord(1),0*elementcoord(2),elementcoord(3),0,SA(element,3)/1e10,0,'r')
%           end
%           
%           for i = 1:length(elelist2)
%                 hold on
%                 element = elelist2(i);
%                 coords = mesh.elementcreate(element).coordinates();
%                 elementcoord = mesh.elementcoord(element);
%                 quiver3(elementcoord(1),200/1e3,elementcoord(3),0,SA(element,3)/1e10,0,'r')
%           end
%           max(SA(elelist,2))
%           for i = 1:length(elelist3)
%                 hold on
%                 element = elelist3(i);
%                 coords = mesh.elementcreate(element).coordinates();
%                 elementcoord = mesh.elementcoord(element);
%                 quiver3(elementcoord(1),100/1e3,elementcoord(3),0,SA(element,3)/1e10,0,'r')
%           end
%           for i = 1:length(elelist4)
%                 hold on
%                 element = elelist4(i);
%                 coords = mesh.elementcreate(element).coordinates();
%                 elementcoord = mesh.elementcoord(element);
%                 quiver3(elementcoord(1),50/1e3,elementcoord(3),0,SA(element,2)/1e10,0,'b')
%           end
%             hold off
%             
%             
%                       %
%           
%           case_NX.plot(case_NX.get('MaxStressIndEle'))
%             dis = case_NX.get('displacements');
%             Dis = dis.xyzout();
%             scale = 10;
%             hold on
%             mesh = case_NX.get('mesh');
%             mesh = MeshClass(mesh.get('coordinates')+Dis*scale,mesh.get('connections'),mesh.get('material'));
% 
% 
%             mesh.plot('g');
%           for i = 1:length(elelist5)
%                 hold on
%                 element = elelist5(i);
%                 coords = mesh.elementcreate(element).coordinates();
%                 elementcoord = mesh.elementcoord(element);
%                 if elementcoord(3)>0.05
%                     aux = abs(SA(element,2)/1e10);
%                 else aux = -abs(SA(element,2)/1e10);
%                 end
%                 quiver3(elementcoord(1),150/1e3,elementcoord(3),0,aux,0,'r')
%           end
%           max(SA(elelist,2))
%           for i = 1:length(elelist4)
%                 hold on
%                 element = elelist4(i);
%                 coords = mesh.elementcreate(element).coordinates();
%                 elementcoord = mesh.elementcoord(element);
%                 if elementcoord(3)>0.05
%                     aux = abs(SA(element,2)/1e10);
%                 else aux = -abs(SA(element,2)/1e10);
%                 end
%                 quiver3(elementcoord(1),50/1e3,elementcoord(3),0,aux,0,'r')
%           end
%           for i = 1:length(elelist4)
%                 hold on
%                 element = elelist4(i);
%                 coords = mesh.elementcreate(element).coordinates();
%                 elementcoord = mesh.elementcoord(element);
%                 quiver3(elementcoord(1),50/1e3,elementcoord(3),0,SA(element,2)/1e10,0,'b')
%           end
%             hold off
%             
%                                  %
%           
%           case_NX.plot(case_NX.get('MaxStressIndEle'))
%             dis = case_NX.get('displacements');
%             Dis = dis.xyzout();
%             scale = 10;
%             hold on
%             mesh = case_NX.get('mesh');
%             mesh = MeshClass(mesh.get('coordinates')+Dis*scale,mesh.get('connections'),mesh.get('material'));
% 
% 
%             mesh.plot('g');
%           for i = 1:length(elelist5)
%                 hold on
%                 element = elelist5(i);
%                 coords = mesh.elementcreate(element).coordinates();
%                 elementcoord = mesh.elementcoord(element);
%                 if elementcoord(3)>0.05
%                     aux = -abs(SA(element,3)/1e10);
%                 else aux = abs(SA(element,3)/1e10);
%                 end
%                 quiver3(elementcoord(1),150/1e3,elementcoord(3),0,aux,0,'r')
%           end
%           max(SA(elelist,2))
%           for i = 1:length(elelist4)
%                 hold on
%                 element = elelist4(i);
%                 coords = mesh.elementcreate(element).coordinates();
%                 elementcoord = mesh.elementcoord(element);
%                 if elementcoord(2)>0.05
%                 if elementcoord(3)>0.05 && elementcoord(2)<0.05
%                     aux = -abs(SA(element,3)/1e10);
%                 else aux = abs(SA(element,3)/1e10);
%                 end
%                 end
%                 quiver3(elementcoord(1),50/1e3,elementcoord(3),0,aux,0,'r')
%           end
%           for i = 1:length(elelist4)
%                 hold on
%                 element = elelist4(i);
%                 coords = mesh.elementcreate(element).coordinates();
%                 elementcoord = mesh.elementcoord(element);
%                 quiver3(elementcoord(1),50/1e3,elementcoord(3),0,SA(element,2)/1e10,0,'b')
%           end
%             hold off
%             
%             
% case_NX.plot(case_NX.get('MaxStressIndEle'))
% 
%             %
% 
%            