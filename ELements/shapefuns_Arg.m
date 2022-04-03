function [B,G,D,Bx,By,Bxx,Byy,Bxy] = shapefuns_Arg(x,y,emap,vmap)

x = x(:); 
y = y(:);
z = zeros(size(x));
o = ones(size(x));

[B,G,D,Bx,By,Bxx,Byy,Bxy] = poly2d_Arg(x,y,z,o);

A = shape2dArg(emap,vmap);

B  = B/A;
G  = G/A;
D  = D/A;
Bx = Bx/A;
By = By/A;
Bxx = Bxx/A;
Byy = Byy/A;
Bxy = Bxy/A;

D  = diag((1./x))*Bx+D;  %cylindrical

end







