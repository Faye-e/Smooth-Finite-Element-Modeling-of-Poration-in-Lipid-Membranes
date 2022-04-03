function  varargout = geval(fun,direct,varargin)

global fe
if ~iscell(direct),   direct = {direct};  end
n = numel(varargin);
V = cell(n,fe.el);
for i = 1:fe.el, for j= 1:n, V{j,i} = varargin{j}(:,i); end; end

f = feval(fun,direct{:},V{:,1});
m = numel(f);
Q = cell(m,fe.el);
Q(:,1) = f;

if fe.par == 0
    for i = 2:fe.el
        Q(:,i) = feval(fun,direct{:},V{:,i});
    end
elseif fe.par == 1
    parfor i = 2:fe.el
        Q(:,i) = feval(fun,direct{:},V{:,i});
    end
end

varargout = cell(1,m);
for k = 1:m, q = [Q{k,:}]; varargout{k}  = q(:); end
end

