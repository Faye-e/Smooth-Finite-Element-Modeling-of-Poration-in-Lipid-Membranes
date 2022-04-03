function [R,U,J]=AllmyRezJac(pa,W,o,l,Iae, r,z,Q,B,Bx,By,Bxx,Byy,Bxy,dtrm)
               
B = B{1};
Br = Bx{1};
Bz = By{1};
Brr = Bxx{1};
Bzz = Byy{1};
Brz = Bxy{1};

%==========================================================================
u   = B*Q;
ur  = Br*Q;
uz  = Bz*Q;
urr = Brr*Q;
uzz = Bzz*Q;
urz = Brz*Q;

%==========================================================================
L = [B;Br;Bz;Brr;Bzz;Brz];

S = kron(speye(6),spd( dtrm*W(:) )); % r should be considered in the generated function

%==========================================================================

[q,qr,qz] = q_fun(pa.r0,pa.a,pa.b,pa.m*sqrt(2),r,z);

[R,Da,J11,J12,J13,J14,J15,J16,J22,J23,J24,J25,J26,J33,J34,J35,J36,J44,J45,J46,J55,J56,J66]=...
    myRezJac(o,l,pa.k,pa.c1,pa.c2,pa.c3,pa.kG,Iae,pa.a0,pa.e,pa.m,pa.l,r,u,ur,uz,urr,uzz,urz,q,qr,qz);

J=[ spd(J11),spd(J12),spd(J13),spd(J14),spd(J15),spd(J16);
    spd(J12),spd(J22),spd(J23),spd(J24),spd(J25),spd(J26);
    spd(J13),spd(J23),spd(J33),spd(J34),spd(J35),spd(J36);
    spd(J14),spd(J24),spd(J34),spd(J44),spd(J45),spd(J46);
    spd(J15),spd(J25),spd(J35),spd(J45),spd(J55),spd(J56);
    spd(J16),spd(J26),spd(J36),spd(J46),spd(J56),spd(J66)];

R = (L'*S*R);
U = (L'*S*Da);
J = (L'*(S*J)*L);  J = (J+J')/2;

%==========================================================================
R  = {R,U,J(:)};
end
function sa = spd(a)
n=numel(a);
I = 1:n;
sa = sparse(I,I,a,n,n);
end