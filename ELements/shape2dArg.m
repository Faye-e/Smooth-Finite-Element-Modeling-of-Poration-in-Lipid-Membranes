function A = shape2dArg(emap,vmap)

filename = 'Afun_Arg';
%--------------------------------------------------------------------------

[a,b,c]=deal(emap{:});
a1=a(1); a2=a(2); a3=a(3); a4=a(4);
b1=b(1); b2=b(2); b3=b(3); b4=b(4);
c1=c(1); c2=c(2); c3=c(3); c4=c(4);

d = vmap(:);
d1=d(1); d2=d(2); d3=d(3); d4=d(4); d5=d(5); d6=d(6);

A = feval(filename,...
    a1,a2,a3,a4,...
    b1,b2,b3,b4,...
    c1,c2,c3,c4,...
    d1,d2,d3,d4,d5,d6);

end
