function gdof = findof(gen,ldof,dof)
% finding global degree of freedom
% syntax:   gdof = findof(gen,ldof,dof)
% gen  = Numbers corresponding to geometry: nodes edges
% ldof = degree of freedom Or generalized coordinate
% dof  = number of dofs related corresponding to nodes or edges, ...
% gdof = global degree of freedom
gdof = bsxfun(@plus, (gen(:)-1)*dof, ldof(:)').'; 
gdof = gdof(:);
end