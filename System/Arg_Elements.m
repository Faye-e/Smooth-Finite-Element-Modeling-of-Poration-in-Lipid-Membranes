%% 
[p,W] = GaussPoints(gg); 
S = diag(W);
[emap,tan,le] = affineMap1D(cord,mnodesXvnodes,mnodesXelems);
[J,invJ,dtrm,b,gc] = affineMap2D(vcord,vnodesXelems,p);
vmap =[J;b];

%% Arg

[B,G,D,Bx,By,Bxx,Byy,Bxy] = geval('Arg',{[]},gc{1},gc{2}, emap,vmap);

m = numel(W);
B2cell = @(x) mat2cell(reshape(x,m,[]),m,(fe.ndof_el + fe.mdof_el)*ones(1,fe.el));

B   = B2cell(B); 
Bx  = B2cell(Bx); 
By  = B2cell(By);
Bxx = B2cell(Bxx);
Byy = B2cell(Byy);
Bxy = B2cell(Bxy);

