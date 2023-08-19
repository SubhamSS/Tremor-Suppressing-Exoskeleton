function mpc_inp = mpc_no_BMFLC(AA,BB,P,Q,R,hor,X_in,F_cons,path2, uMax, uMin, xx_lb, xx_ub, hh)     
% equality constraints
    %ny1 = size(AA,2); 
    nu =size(BB,2); ny = size(AA,1);
    bigA = zeros(ny*hor,ny*hor);
    for ii = 1:ny*hor
        bigA(ii,ii)=1;
    end
    for ii = 1:hor-1
        for jj = 1:ny
            for kk = 1:ny
                bigA(ny*ii+jj,ny*(ii-1)+kk)=-AA(jj,kk); 
            end
        end
    end
%     bigA = [AA, zeros(ny,ny1*(hor-1))];
%     for ii = 1:hor-1
%         bigA = [bigA ; zeros(ny,ny1*(ii-1)),-AA,AA,zeros(ny,ny1*(hor-ii-1))];
%     end
%     
    bigB = zeros(ny*hor,nu*hor);
    for ii = 0:hor-1
        for jj = 1:ny
            for kk = 1:nu
                bigB(ny*ii+jj,nu*ii+kk)=-BB(jj,kk); 
            end
        end
    end
    Eq_RHS = [AA*X_in;zeros(ny*(hor-1),1)] + repmat(F_cons,hor,1);
    bigEq = [bigA,bigB];
    
    
    % Cost function Matrix 
    Qbar = zeros(ny*hor,ny*hor);
    for ii = 0:hor-1
        for jj = 1:ny
            for kk = 1:ny
                if ii<hor-1
                    Qbar(ny*ii+jj,ny*ii+kk)=Q(jj,kk); 
                else
                    Qbar(ny*ii+jj,ny*ii+kk)=P(jj,kk); 
                end
            end
        end
    end
    
    Rbar = zeros(nu*hor,nu*hor);
    for ii = 0:hor-1
        for jj = 1:nu
            for kk = 1:nu
                    Rbar(nu*ii+jj,nu*ii+kk)=R(jj,kk); 
            end
        end
    end
    QR_bar = [Qbar,zeros(ny*hor,nu*hor);zeros(nu*hor,ny*hor),Rbar];
    
    % Inequality constraints
    
    u_lb = repmat(uMin,hor,1);
    u_ub = repmat(uMax,hor,1);

    %xx_lb = [-1000;-1000;-1000;-1000]; %taken randomly for now
    x_lb = repmat(xx_lb,hor,1);
    %xx_ub = [10;10;1000;1000];
    x_ub = repmat(xx_ub,hor,1);
   
    Aineq = [eye(hor*(ny+nu));-eye(hor*(ny+nu))];
    Bineq = [x_ub;u_ub;-x_lb;-u_lb];

    %qp solve
    ff = [];
    for ii = 1:1:hor
        ff = [ff,-1*path2(:,ii)'*Q];
    end
    ff = [ff, zeros(1,nu*hor)];

    ff1 = ff';
    [x1_u1,~,~]=qpSWIFT(sparse(QR_bar),ff1,sparse(bigEq),Eq_RHS,sparse(Aineq),Bineq);


    mpc_inp = x1_u1(4*hor+1:4*hor+20);
   
end

