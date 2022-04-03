function [I,J,L] = gindicesComplex(contab,dof)
if ~iscell(contab)
    
    t = contab;
    n = size(t,1)*dof;
    
    t = bsxfun(@plus, (t(:).'-1)*dof, (1:dof).');
    L = t(:);     %assembling vectors
    
    I = repmat( reshape(t,n,[]), [n,1]);
    I = I(:);
    
    J = ( t(:)*ones(1,n) ).';
    J = J(:);
    
else    
    if ~iscell(dof), dof = {dof(1),dof(2)}; end    
    for j = 1:2
        ct = contab{j};
        if ~iscell(ct), ct = {ct}; end
        d  = dof{j};
        ta = [];
        s = 0;
        for i = 1:numel(ct)
            t = ct{i};
            n = size(t,1)*d(i);
            t = t(:);                      
            t = bsxfun(@plus, (t'-1)*d(i), (1:d(i)).');            
            ta = [ta; s+reshape(t,n,[])];
            s = max(max(ta)); 
        end
        tt{j} = ta;
    end
    
    L = tt{1}(:);     %assembling vectors
    
    I = repmat(tt{1},[size(tt{2},1),1]);
    I = I(:);
    
    J = ( tt{2}(:)*ones(1,size(tt{1},1)) ).';
    J = J(:);
end

end


