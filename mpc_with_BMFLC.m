function mpc_inp = mpc_with_BMFLC(AA,BB,P,Q,R,S_gain,hor,X_in,F_cons,path2, uMax, uMin, xx_lb, xx_ub, tr_jac,hh)     
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
    Eq_RHS2 = [Eq_RHS;zeros(64*1,1)]; %$%
    bigEq = [bigA,bigB];
    B_times_f = BB*tr_jac;
    bigbigEq = [bigEq,zeros(ny*hor,64*1);zeros(64*1,size(bigEq,2)),zeros(64*1,64*1)];
    for i=1:1:hor
        bigbigEq((i-1)*ny+1:(i-1)*ny+ny,hor*(ny+nu)+1:hor*(ny+nu)+64)=B_times_f;
    end
    
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
    S_bar = S_gain*eye(64*1);
    QRS_bar = [QR_bar,zeros(nu*hor+ny*hor,64*1);...
        zeros(64*1,nu*hor+ny*hor),S_bar];
    
    % Inequality constraints
    
    u_lb = repmat(uMin,hor,1);
    u_ub = repmat(uMax,hor,1);

    %xx_lb = [-1000;-1000;-1000;-1000]; %taken randomly for now
    x_lb = repmat(xx_lb,hor,1);
    %xx_ub = [10;10;1000;1000];
    x_ub = repmat(xx_ub,hor,1);
   
    Aineq = [eye(hor*(ny+nu));-eye(hor*(ny+nu))];
    Aineq2 = zeros(hor*(ny+nu),1*64);
    for i = 1:1:1
        Aineq2(hor*ny+(i-1)*nu+1:hor*ny+(i-1)*nu+2,1:64)=tr_jac;
    end
    Aineq3 = zeros(hor*(ny+nu),1*64);
    for i = 1:1:1
        Aineq3(hor*ny+(i-1)*nu+1:hor*ny+(i-1)*nu+2,1:64)=-tr_jac;
    end
    Aineqr = [Aineq2;Aineq3];
    Aineq_fin = [Aineq,Aineqr];
    Bineq = [x_ub;u_ub;-x_lb;-u_lb];

    %qp solve
    ff = [];
    for ii = 1:1:hor
        ff = [ff,-1*path2(:,ii)'*Q];
    end
    ff = [ff, zeros(1,nu*hor), zeros(1,64*1)];

    ff1 = ff';
    [x1_u1,~,~]=qpSWIFT(sparse(QRS_bar),ff1,sparse(bigbigEq),Eq_RHS2,sparse(Aineq_fin),Bineq);

    mpc_inp = zeros(20,1);
    %x1_u1(ny*hor+1:ny*hor+20);
    inp2 = zeros(20,1);
    for i = 1:1:10
           inp2(2*i-1:2*i,1)=tr_jac*x1_u1((ny+nu)*hor+1:(ny+nu)*hor+64);
    end
    for i = 1:1:10
           mpc_inp(2*i-1:2*i,1)=x1_u1(ny*hor+(i-1)*nu+1:ny*hor+(i-1)*nu+2);
    end
    
end

