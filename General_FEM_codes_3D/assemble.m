function varargout = assemble(locals,contab,dof)

% [I,J,L] = gindices(contab,dof);
[I,J,L] = gindicesComplex(contab,dof);

if ~iscell(locals)
    if numel(L) == numel(locals);     varargout = {sparse(L,1,locals)};
    else                              varargout = {sparse(I,J,locals)};
    end    
else
    m = numel(locals);                varargout = cell(1,m);
    for i = 1:m
        if numel(L)==numel(locals{i});varargout{i} = sparse(L,1,locals{i});
        else                          varargout{i} = sparse(I,J,locals{i});
        end
    end
end

end
