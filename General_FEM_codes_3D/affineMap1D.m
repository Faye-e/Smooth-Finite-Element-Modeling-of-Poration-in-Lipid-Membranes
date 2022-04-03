function [emap,tan,le] = affineMap1D(cord,edges_vnodes,edges_elems)

ne = size(edges_vnodes,1); % number of edges
emap = cell(ne,1);  %from face 1D coordinate s to global coordinates
le = zeros(ne,1);
tan = cell(ne,1);

for i = 1:ne
    vce = cord(edges_vnodes(i,:),:)';
    Ab = [vce(:,2)-vce(:,1), vce(:,1)];   %[A,b]
    tan0 = Ab(:,1);
    l = norm(tan0);
      
    le  (i) = l;
    emap{i} = Ab;
    tan {i} = tan0/l;
end
%--------------------------------------------------------------------------
% Allocating data to elements
le   = le  (edges_elems); 
emap = emap(edges_elems);
tan  = tan (edges_elems);

end
