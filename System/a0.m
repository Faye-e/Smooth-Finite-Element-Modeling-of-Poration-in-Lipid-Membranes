function Iae = a0(Q,vnodesXelems,mnodesXelems,r,z,B,Bx,By,W,dtrm)

global fe pa 
 
Qn = Q(gdofs(vnodesXelems,fe.ndof)); 
Qm = Q(gdofs(mnodesXelems,fe.mdof) + fe.ngdof);
Q = [Qn;Qm]; 


Iae = geval('AllmyAeVe',{pa,W},r,z,Q,B,Bx,By,dtrm);
Iae = sum(Iae); % ON ENTIRE DOMAIN