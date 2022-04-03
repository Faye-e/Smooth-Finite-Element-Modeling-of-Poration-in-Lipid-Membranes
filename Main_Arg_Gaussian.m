clc; clear; warning on; 
path(pathdef);
addpath('General_FEM_codes_3D','Solvers','System','ELements', 'VTKformat');
global fe pa par1 par2
 

%% Parametersl
pa.c0 = 0;
pa.a = 0;
aui = 0;
initial = 'sphere';

gg  = 13;
h0  = 0.03;
pa.e = 0.1;
pa.m = pa.e;


par2 = linspace(1,1,100);


par1 = 0.2*ones(size(par2)); 
par  = 1:numel(par2);

pa.r0 = 0.3;
pa.b  = par1(1);

rui = 0.3; 
bui = -0.1; 

pa.k  = 1;
pa.kG = -0.1*0;
pa.l = par2(1);              

pa.c1 = 1e3;
pa.c2 = 1e3;
pa.c3 = 1e3;
pa.c4 = 1e3;
pa.c5 = 1e3;

d = [0   -0.8;
     0.8  0.8];

%% Generating the Mesh

pa.wR = prod(d(2,:)-d(1,:));
pa.aR = d(2,1)-d(1,1); 

g = decsg([2 4  d(1) d(2) d(2) d(1)  d(3) d(3) d(4) d(4) ]');
[p,e,t] = initmesh(g,'Hmax',h0);
for i = 1:0, [p,e,t] = refinemesh(g,p,e,t,'regular'); end
p = jigglemesh(p,e,t,'Opt','minimum');
vcord = p'; vnodesXelems = t(1:3,:);

%%
MeshAndBoundaries_Arg

%% Initial guess               

switch initial
    case 'sphere' 
        g_u = @(r,z) sqrt((r-aui).^2+(z-bui).^2);
end

u_ind = findof(1:fe.n, 1, fe.ndof);
Q = sparse(fe.tdof,1);
Q(u_ind) = tanh((rui - g_u(vR,vZ))/pa.e/sqrt(2));
Q(bound) = bval;

%% plot initial guess
un0 = full(Q(u_ind));

q = q_fun(pa.r0,pa.a,pa.b,pa.m*sqrt(2),full(vR),full(vZ));

figure(5); clf;  

subplot(1,3,1)
trisurf(vnodesXelems',full(vR),full(vZ),un0); hold on;
trisurf(vnodesXelems',full(vR),full(vZ),q);  
title([ 'Guess - #Elements = ' num2str(size(vnodesXelems,2)) ]);
view([0 90]);
axis equal

subplot(1,3,2)
tr =  0.1;
un1 = un0;
Id = (q > tr) | un1<-tr | un1>tr;
un1( Id) = 0;
un1(~Id) = 1;

lvr = vR(~Id);
lvz = vZ(~Id);
trisurf(vnodesXelems',full(vR),full(vZ),un1); hold on;
shading interp
title([ 'Guess - #Elements = ' num2str(size(vnodesXelems,2)) ]);
view([0 90]);
axis equal

subplot(1,3,3)
trisurf(vnodesXelems',full(vR),full(vZ),un0); hold on;
trisurf(vnodesXelems',full(vR),full(vZ),q);
shading interp

rr = sqrt(lvr.^2+lvz.^2);
tt = atan(lvz./lvr);
tab = sortrows([rr(:),tt(:)],2);
rr = tab(:,1);
tt = tab(:,2);
p = polyfit(tt,rr,7);
tt = linspace(min(tt),max(tt),100);
rr = polyval(p,tt);
lvr = rr.*cos(tt);
lvz = rr.*sin(tt);
plot3(lvr,lvz,ones(size(lvz)),'m-','linewidth',4)


view([0 90]);
title([ 'Guess - #Elements = ' num2str(size(vnodesXelems,2)) ]);
axis equal

pause
%%
Arg_Elements

%%

pa.a0 = a0(Q,vnodesXelems,mnodesXelems,gc{1},gc{2},B,Bx,By,W,dtrm); 

%% Nonlinear solution

fun='NonlinSys_shmo';

opt = arcset('newitr',300,'newer',0.05);
%                
Ud = newcont(fun,Q(domain),par,opt,15,...
    domain,bound,bval,vnodesXelems,mnodesXelems,W,gc{1},gc{2},...
    B,Bx,By,Bxx,Byy,Bxy,dtrm,...
    domain,bound,bval,vnodesXelems,vR,vZ);

for ii=1:size(Ud,2)
    U = sparse(domain,1,Ud(:,ii),fe.tdof,1);
    U(bound) = bval;
    UU(:,ii) = U;
end


