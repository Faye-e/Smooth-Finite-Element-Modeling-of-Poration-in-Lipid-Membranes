%% coordinates
R = vcord(:,1); 
Z = vcord(:,2);
vR = R; 
vZ = Z;

%% Mesh

% figure(1); clf; hold on
% colormap('parula')
% trisurf(vnodesXelems',full(R),full(Z),zeros(size(R)),'facecolor',[0.5 0.5 0.5]); view([0 90]);
% title([ '#Elements = ' num2str(size(vnodesXelems,2)) ]);
% axis equal; axis tight
% drawnow

%% Connectivity Tables
oc = orderCode2D;
    
[mnodesXelems, mnodesXvnodes] = partsFinder(oc.e,vnodesXelems');
[bmnodes, bnodes] = getbounds3D(mnodesXelems, mnodesXvnodes);

mcord = [mean(vR(mnodesXvnodes),2),mean(vZ(mnodesXvnodes),2)];
mR = mcord(:,1);
mZ = mcord(:,2);

cord = [vcord;mcord];
R = cord(:,1); 
Z = cord(:,2);
nodesXelems = [vnodesXelems; mnodesXelems; vnodesXelems];

%%
tg = [vR(mnodesXvnodes)*[-1;1],vZ(mnodesXvnodes)*[-1;1]];  %tangents
Le = sqrt(tg(:,1).^2+tg(:,2).^2);
% tg = diag(1./Le)*tg;              %normalizing
rr = [0 1;-1 0];  %rotation
%rr = eye(2);
nor = (tg*rr')';

%% Number of nodes, Elements, and DOFs
fe = feset2d_Arg(vcord,mcord,vnodesXelems,mnodesXelems);
fe.shmo = 1;
fe.par = 1;

if fe.par == 1 && isempty(gcp('nocreate')), parpool('local'); end
disp(fe.tdof)

%% Findig Boundaries
nn = unique(vnodesXelems);             %node numbers (condiders the index of common nodes only one time)
er = 1e-6;
clear B
B.r = nn( abs(vR(nn)-max(vR(nn)))<er);
B.l = nn( abs(vR(nn)-min(vR(nn)))<er);
B.t = nn( abs(vZ(nn)-max(vZ(nn)))<er);
B.b = nn( abs(vZ(nn)-min(vZ(nn)))<er);

nn = unique(mnodesXelems);
Bm.r = nn( abs(mR(nn)-max(mR(nn)))<er);
Bm.l = nn( abs(mR(nn)-min(mR(nn)))<er);
Bm.t = nn( abs(mZ(nn)-max(mZ(nn)))<er);
Bm.b = nn( abs(mZ(nn)-min(mZ(nn)))<er);

%% Plot the boundaries of the mesh
% figure(1); hold on;
% plot(vR(B.l),vZ(B.l),'o w','markerfacecolor','b','markersize',10);
% plot(vR(B.r),vZ(B.r),'o w','markerfacecolor','b','markersize',10);
% plot(vR(B.t),vZ(B.t),'o w','markerfacecolor','r','markersize',10);
% plot(vR(B.b),vZ(B.b),'o w','markerfacecolor','r','markersize',10);
% %==========================================================================
% plot(mR(Bm.l),mZ(Bm.l),'s w','markerfacecolor','b','markersize',10);
% plot(mR(Bm.r),mZ(Bm.r),'s w','markerfacecolor','b','markersize',10);
% plot(mR(Bm.t),mZ(Bm.t),'s w','markerfacecolor','r','markersize',10);
% plot(mR(Bm.b),mZ(Bm.b),'s w','markerfacecolor','r','markersize',10);
% xlabel('R'); ylabel('Z'); 
% axis equal

%% Findig Boundaries

er = 1e-6;
mm = [max(vcord);min(vcord)];

bc.r = bnodes( abs(R(bnodes)-mm(1,1))<er);
bc.l = bnodes( abs(R(bnodes)-mm(2,1))<er);
bc.t = bnodes( abs(Z(bnodes)-mm(1,2))<er);
bc.b = bnodes( abs(Z(bnodes)-mm(2,2))<er);

ibm = unique(bmnodes);

%% u ux uy uxx uyy uxy

bound = [
    findof( bc.r , [1 2 3 5 6], fe.ndof);   % u ux uy uyy uxy
    findof( bc.l , [  2     6], fe.ndof);   %   ux        uxy
    findof( bc.t , [1 2 3 4 6], fe.ndof);   % u ux uy uxx uxy
    findof( bc.b , [1 2 3 4 6], fe.ndof);   % u ux uy uxx uxy
    %-------------------------------------------------------------
    findof( ibm ,1, fe.mdof) + fe.ngdof;   % normal at mid point
    %-------------------------------------------------------------
    ];

domain = setdiff(1:fe.tdof,bound)';

%%

Rr = R(bc.r);  Zr = Z(bc.r);
Rl = R(bc.l);  Zl = Z(bc.l);
Rt = R(bc.t);  Zt = Z(bc.t);
Rb = R(bc.b);  Zb = Z(bc.b);

u   = @(r,z) -ones(numel(r),1);
ur  = @(r,z) zeros(numel(r),1);
uz  = @(r,z) zeros(numel(r),1);
urr = @(r,z) zeros(numel(r),1);
uzz = @(r,z) zeros(numel(r),1);
urz = @(r,z) zeros(numel(r),1);

uvalr = [u(Rr,Zr) ur(Rr,Zr) uz(Rr,Zr) uzz(Rr,Zr) urz(Rr,Zr)]';
uvall = [         ur(Rl,Zl)                      urz(Rl,Zl)]';
uvalt = [u(Rt,Zt) ur(Rt,Zt) uz(Rt,Zt) urr(Rt,Zt) urz(Rt,Zt)]';
uvalb = [u(Rb,Zb) ur(Rb,Zb) uz(Rb,Zb) urr(Rb,Zb) urz(Rb,Zb)]';
%------------------------------
M = [];
M(1,:) = ur( mR(ibm), mZ(ibm) )';
M(2,:) = uz( mR(ibm), mZ(ibm) )';
DOFm = sum(M.*nor(:,ibm));

bval = [
    uvalr(:);
    uvall(:);
    uvalt(:);
    uvalb(:);
    %------------------------------
    DOFm'
    ];
