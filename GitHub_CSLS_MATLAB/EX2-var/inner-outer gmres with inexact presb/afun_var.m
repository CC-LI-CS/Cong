%% Find afun=YAx        (minres(A,b,tol,maxit,gfun))

function vout = afun_var(alpha1,beta1,Cxpmat,Cypmat,Cdmat,dim,M,Ms,N,vin)            

%% 
x1=vin(1:Ms);
x2=vin(Ms+1:end);

%% La*x1
vin1=reshape(x1,[M,M,N]);                                    

vout1=reshape(repmat(Cdmat,[1,1,N]).*vin1...
    +repmat([-Cxpmat;zeros(1,M)],[1,1,N]).*[vin1(2:end,:,:);zeros(1,M,N)]...
    +repmat([zeros(1,M);-Cxpmat],[1,1,N]).*[zeros(1,M,N);vin1(1:end-1,:,:)]...
    +repmat([-Cypmat,zeros(M,1)],[1,1,N]).*[vin1(:,2:end,:),zeros(M,1,N)]...
    +repmat([zeros(M,1),-Cypmat],[1,1,N]).*[zeros(M,1,N),vin1(:,1:end-1,:)],M^dim,N);

vecout=vout1(:,1:N);
vout1=vecout(:);

%% La*x2
vin2=reshape(x2,[M,M,N]);                                    
vout2=reshape(repmat(Cdmat,[1,1,N]).*vin2...
    +repmat([-Cxpmat;zeros(1,M)],[1,1,N]).*[vin2(2:end,:,:);zeros(1,M,N)]...
    +repmat([zeros(1,M);-Cxpmat],[1,1,N]).*[zeros(1,M,N);vin2(1:end-1,:,:)]...
    +repmat([-Cypmat,zeros(M,1)],[1,1,N]).*[vin2(:,2:end,:),zeros(M,1,N)]...
    +repmat([zeros(M,1),-Cypmat],[1,1,N]).*[zeros(M,1,N),vin2(:,1:end-1,:)],M^dim,N);

vecout=vout2(:,1:N);
vout2=vecout(:);

%% 
vout=[beta1*x1+alpha1*x2+vout2;alpha1*x1-beta1*x2+vout1];
% aa=size(vout)
end