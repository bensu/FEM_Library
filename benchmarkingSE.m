clear all
dim = 2;
coordnum = 2;
plot_on = true;
%% Mesh & Material
a = 2; m = dim;
b = 1; n = dim;
c = 3; p = dim;

E = 2;                                  % Youngs modulus
nu = 0;                                 % Poisson's ratio
rho = 7840;                           % Density of the plate
sides = [a,b,c];

%% Fist Mesh
mesh = MeshClass.meshgen(sides,[m n p],E,nu,rho);

bcface = Face.new_face(mesh,coordnum,mesh.mincoordval(coordnum));
middleface = Face.new_face(mesh,coordnum,mesh.maxcoordval(coordnum));
mesh.set('nodesout',[bcface.get('nodelist');middleface.get('nodelist')]);


%% Second Mesh
mesh2 = MeshClass.meshgen(sides,[m n p],E,nu,rho);

bcface = Face.new_face(mesh2,coordnum,mesh2.mincoordval(coordnum));
middleface = Face.new_face(mesh2,coordnum,mesh2.maxcoordval(coordnum));
mesh2.set('nodesout',[bcface.get('nodelist');middleface.get('nodelist')]);

shift_vector = zeros(1,3);
shift_vector(coordnum) = mesh.maxcoordval(coordnum);
mesh2.shift(shift_vector);

mesh.get('coordinates')
mesh2.get('coordinates')

%% FEM CASE

patch = FemCase.case_from_SElist([mesh mesh2]);
patch.get('mesh').coordinates
patch.get('mesh').plot('g')

%% BC
bc_face = Face.new_face(patch.get('mesh'),coordnum,patch.get('mesh').mincoordval(coordnum));
patch.get('bc').simplysupportedface(bc_face,coordnum,patch.get('mesh'));

%% Loads
coordload = coordnum;
q = 2;
vector = zeros(3,1);
vector(coordload) = q;
load_face = Face.new_face(patch.get('mesh'),coordnum,patch.get('mesh').maxcoordval(coordnum));
patch.get('loads').qinface(patch.get('mesh'),load_face,vector);


%% CASE
patch.solve_condensed()
patch.getMaxStressEle()
%subplot(3,1,coordnum)
if plot_on
    patch.plot()
end

%% check

dis = patch.get('displacements');
maxval = dis.maxvalue();
tol = 1e-12;
assertEqual(abs(maxval-(q*sides(coordnum)/E))<tol,true)