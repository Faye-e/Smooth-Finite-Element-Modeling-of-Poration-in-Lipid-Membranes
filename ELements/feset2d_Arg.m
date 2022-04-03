function fe = feset2d_Arg(vcord,mcord,vnodes_elems,mnodes_elems)
% Tables of connectivities  [dofs  X elems]
% coordinates               [nodes X cords]
%--------------------------------------------------------------------------
fe.el   = size(vnodes_elems,2);   % total number of elements
fe.n    = size(vcord,1);          % total number of number of vnodes
fe.m    = size(mcord,1);          % total number of number of mnodes
fe.tn   = fe.n + fe.m;
fe.n_el = size(vnodes_elems,1);    % number of vnodes in each element
fe.m_el = size(mnodes_elems,1);    % number of mnodes in each element
%--------------------------------------------------------------------------
fe.ndof = 6;        % number of dofs for each vnode 
fe.mdof = 1;        % number of dofs for each mnode
%--------------------------------------------------------------------------
fe.ngdof = fe.n*fe.ndof;   
fe.mgdof = fe.m*fe.mdof;  
fe.tdof  = fe.ngdof + fe.mgdof;     % TOTAL dofs
%---------------------------------
fe.ndof_el = fe.n_el*fe.ndof;
fe.mdof_el = fe.m_el*fe.mdof;
fe.tdof_el = fe.ndof_el + fe.mdof_el;  % total dofs in each element
end
