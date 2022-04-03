function [bfaces,bnodes] = getbounds3D(faces_elems, faces_nodes)

faces_rep = faces_elems(:);

bfaces =  find( sum(bsxfun(@eq,unique(faces_rep ),faces_rep'),2)==1 );

bfaces_nodes = faces_nodes(bfaces,:);
bnodes = unique(bfaces_nodes(:));

end
