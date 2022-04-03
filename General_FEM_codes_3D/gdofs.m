function gd = gdofs(contab,dof)
    n = size(contab,1)*dof;
    t = bsxfun(@plus, (contab(:).'-1)*dof, (1:dof).');
    gd =  reshape(t,n,[]);    
end