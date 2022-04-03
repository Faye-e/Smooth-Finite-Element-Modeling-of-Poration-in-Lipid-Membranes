function out = AllmyAeVe(pa,W,  r,z,Q,B,Br,Bz,dtrm)
B = B{1};                 
Br = Br{1};
Bz = Bz{1};

W = dtrm*W(:)';

u  = B*Q;
ur = Br*Q;
uz = Bz*Q;

q = q_fun(pa.r0,pa.a,pa.b,pa.m*sqrt(2),r,z);

ae = myAeVe(pa.e,pa.m,r,u,ur,uz,q);
ae = W*ae;

out  = {ae};

end

