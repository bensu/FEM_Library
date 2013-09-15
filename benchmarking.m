clear all
dim = 2;
coordnum = 1;
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

%% BC
valuebc = 0;
face = Face.new_face(mesh,coordnum,valuebc);
bc = BC(mesh.sdof());
bc.simplysupportedface(face,coordnum,mesh);

%% LOADS
coordload = coordnum;

q = 2;
valueload = sides(coordload);
vector = zeros(3,1);
vector(coordload) = q;
face = Face.new_face(mesh,coordload,valueload);
loads = Loads(mesh.sdof());
loads.qinface(mesh,face,vector);

%% CASE
patch = FemCase(mesh,mesh.sdof(),mesh.sdof());
patch.set('loads',loads);
patch.set('bc',bc);
patch.solve();
patch.getMaxStressEle()
SA = patch.get('StressArrayEle');
%subplot(3,1,coordnum)
if plot_on
    patch.plot()
end

%% check

dis = patch.get('displacements');
maxval = dis.maxvalue();
tol = 1e-12;
assertEqual(abs(maxval-(q*sides(coordnum)/E))<tol,true)