function [parts_elems, parts_vnodes] = partsFinder(ocp,elems_vnodes)
el = size(elems_vnodes,1);
[n,m] = size(ocp);

parts_vnodes_rep = zeros(el*n,m); 
for i =1:n
    parts_vnodes_rep(i:n:end,:) = elems_vnodes(:,ocp(i,:));   
end
 
[~,ind_uni,ind_rep] = unique( sort(parts_vnodes_rep,2), 'rows');                                                              

parts_vnodes = parts_vnodes_rep(ind_uni,:); 
parts_elems = reshape(ind_rep, n,[]); 
end