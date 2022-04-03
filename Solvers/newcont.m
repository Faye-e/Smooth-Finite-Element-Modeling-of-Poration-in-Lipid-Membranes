function UU = newcont(fun,X,par,opt,mm,varargin)
global fe pa param


domain       = varargin{mm+1};
bound        = varargin{mm+2};
bval         = varargin{mm+3};
vnodesXelems = varargin{mm+4};
vR           = full(varargin{mm+5});
vZ           = full(varargin{mm+6});

q = q_fun(pa.r0,pa.a,pa.b,pa.m*sqrt(2),vR,vZ);


m = numel(par);
n = size(X,1);
UU = zeros(n,m);
for k = 1:m
    parn = par(k);
        for i = 1:opt.newitr                                %Update by Newton
            warning on
            
            if fe.shmo == 1
                [Fx,U,Jac] = feval(fun,[X;parn],varargin{1:mm});   
                dX = shmo2(Jac,U,Fx);
                disp(['itr = ' num2str(i) ', nomr(dX) = ' num2str(norm(dX))])
                disp(['itr = ' num2str(i) ', nomr(Fx) = ' num2str(norm(Fx))])
                if norm(Fx) <= opt.newer, itr = i-1; break; end
                
            elseif fe.shmo == 0
                [Fx,Jac] = feval(fun,[X;parn],varargin{1:mm});
                disp(['itr = ' num2str(i) ', nomr(F) = ' num2str(norm(Fx))])
                if norm(Fx) <= opt.newer, itr = i-1; break; end
                dX = Jac\Fx;
            end       
             
           X = X-dX;    
           
        %==================================================================       
        U = sparse(domain,1,X,fe.tdof,1);
        U(bound) = bval; 

        un = U(findof(1:fe.n,1,fe.ndof));        

        figure(31); clf
        subplot(1,2,1)
        trisurf(vnodesXelems',vR,vZ,full(un));  hold on      
        trisurf(vnodesXelems',vR,vZ,q);         
        axis equal
        view([0 90])
        shading interp
        title([param ' = ' num2str(parn)])
        
        subplot(1,2,2)
        tr = .1;
        u1 = full(un);
        Id = (q > tr) | u1<-tr | u1>tr;
        u1( Id) = 0;
        u1(~Id) = 1;
        trisurf(vnodesXelems',vR,vZ,u1);
        shading interp
        axis equal
        view([0 90])
        title([param ' = ' num2str(parn)])

        drawnow
        %=================================================================      
        
        end
        
        UU(:,k) = X;
    
        if i ~= opt.newitr
            fprintf('point %g converged in %g iterations \n- param: %g\n',...
                k,i-1,parn);
        else
            beep
            fprintf('point %g did not converge! \n- iteration: %g\n- param: %g\n- ',...
                 k,i-1,parn);
        end
        %==================================================================

end

end
% function X = shmo1(A,U,V,B)
% % Gives X = (A+U*U'+V*V')\B if A=A'.
% % Based on Sherman–Morrison formula.
% % This is faster if A is sparse.
% 
% C=A\[B,V,U];
% M = C(:,2);
% L = C(:,3)-(M'*U)/(1+V'*M)*M;
% X = C(:,1)-(M'*B)/(1+V'*M)*M-(L'*B)/(1+U'*L)*L;
% end

function X = shmo2(A,U,B)
% Gives X = (A+U*U')\B if A=A'.
% Based on Sherman–Morrison formula.
% This is faster if A is sparse.

C=A\[B,U];

X = C(:,1)- (U'*C(:,1))/(1 + U'*C(:,2))*C(:,2);
end
