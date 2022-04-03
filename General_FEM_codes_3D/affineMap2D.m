function [A,iA,detA,b,gpc] = affineMap2D(cord,nodes_elems,p)

% Results in
% T = [A][xi;eta;zeta]+[b]
%-----------------------------------------
R = cord(:,1);  % X coordinates of all nodes in a vector
Z = cord(:,2);
R = R(nodes_elems); % X coordinates of each element (size:3xN)
Z = Z(nodes_elems);

%Trasformation jacobian  
A =[R(2,:)-R(1,:);     
    Z(2,:)-Z(1,:);
    %-------------
    R(3,:)-R(1,:);     
    Z(3,:)-Z(1,:)]; %[A11;A21;A12;A22] %standard vectorizing (:)

b = [R(1,:); Z(1,:)];

%-----------------------------------------
detA = A(1,:).*A(4,:)-A(3,:).*A(2,:);
 
iA = [A(4,:);
     -A(2,:);
     -A(3,:);
      A(1,:)]; %standard vectorizing (:). 

iA = bsxfun(@rdivide,iA,detA); %The inverse of A

if nargin == 3
    % guass points cordinates in each elements
    m = size(p,1); O = ones(m,1);
    gpc{1} = p(:,1)*A(1,:) +p(:,2)*A(3,:) +O*b(1,:); % The mapped X coordinate
    gpc{2} = p(:,1)*A(2,:) +p(:,2)*A(4,:) +O*b(2,:); % The mapped Y coordinate
end
end

% for i = 1:fe.el, Temp(:,i) = diag((reshape(iJ(:,i),[3 3])* reshape(J(:,i),[3 3]))); end
% figure(2);plot(Temp)
