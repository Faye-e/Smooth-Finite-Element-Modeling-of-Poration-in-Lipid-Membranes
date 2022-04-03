function  [R,U,J] = NonlinSys_shmo(Q,domain,bound,bval,vnodesXelems,mnodesXelems,W,r,z,B,Bx,By,Bxx,Byy,Bxy,dtrm)                                     
global fe pa  par1 par2
i = Q(end);

pa.b = par1(i);
pa.l = par2(i);

Q = Q(1:end-1);

%% imposing boundary conditions into initial guess (K linear is not calculated seprately)

Q = sparse(domain,1,Q,fe.tdof,1);
Q(bound) = bval;

%% 
Qn = Q(gdofs(vnodesXelems,fe.ndof)); 
Qm = Q(gdofs(mnodesXelems,fe.mdof) + fe.ngdof);
Q = [Qn;Qm]; 

%% Calculating ae and its drivatives 
Iae = geval('AllmyAeVe',{pa,W},r,z,Q,B,Bx,By,dtrm);
Iae = sum(Iae); % ON ENTIRE DOMAIN
o = zeros(numel(W),1);
l = ones(numel(W),1);

%% Residual and its Jacobian
conectivity = {{vnodesXelems,mnodesXelems}, {vnodesXelems,mnodesXelems} };
DOF = {[fe.ndof fe.mdof],[fe.ndof fe.mdof]};

[R,U,J] = geval('AllmyRezJac',{pa,W,o,l,Iae},r,z,Q,B,Bx,By,Bxx,Byy,Bxy,dtrm);

[R,U,J] = assemble({R,U,J},conectivity,DOF); 
% J is the Jacobian when (Iae) and (Iv) are considered to be constant

%% Imposing boundary condition into res and Jac (elimination approach)

R = R(domain);
U = U(domain);
J = J(domain,domain);
J = (J+J')/2;           % telling matlab that Kt is symmetric 

end
