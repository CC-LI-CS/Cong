function [x]=mgm_inv(b,eigB,y,u_e)

set(0, 'defaultaxesfontsize',20,'defaultaxeslinewidth',1,...
    'defaultlinelinewidth',2,'defaultpatchlinewidth',1.5,...
    'defaulttextfontsize',20,'defaulttextInterpreter','latex');
global LEVEL0;

T=1; mu1=1;mu2=0;  LEVEL0=3; %3
smlist={'Jacobi','Jacobi-C','Vanka','Mass','Mass-shift','Vanka-e','Vanka-e2','GS','RB-GS'};
% lsty ={'x-','d-','^-','v-','>-','<-','p-','h-','.-'};
for kk=[1] %:length(smlist)
    smtype=smlist{kk}; %select compared smoothers
    fprintf('--------------Smoother type: [%s]-------------------\n',smtype)
    for n=y:y
        nt=1;
        %uncomment to select different examples
        %B from Direct PinT solver (heat equation case)
%         alpha1=100;beta1=100;
%         exname='Heat';  dt=T/nt;  %Example 
%         iterlist=[]; ratelist=[];
        for jj=1:nt
            lamj=eigB(jj);
            fineL=n;
            xmin=0;xmax=1;
            ymin=0;ymax=1;
            %Setup operators
            mg=[];
            %define the matrix at each level
            %             nn=2^n; h=(xmax-xmin)/nn;m=nn-1;
            N=2^fineL;
            S0=u_e;
            for Level=fineL:-1:LEVEL0
                nn=2^Level; h=(xmax-xmin)/nn;m=nn-1;
                M=m; % space
               
                xgrid=(xmin:h:(xmax-h)).';
                xpgrid=xgrid+0.5*h;
                xgrid=xgrid(2:end);
                ygrid=(ymin:h:(ymax-h)).';
                ypgrid=ygrid+0.5*h;
               
                ygrid=ygrid(2:end);
                
                %% bar(beta)
                scalx=1/h^2; scaly=1/h^2;
                [Cxpmat,Cypmat,Cdmat]=getcoefmat(scalx,scaly,xgrid,xpgrid,ygrid,ypgrid);%get stepsize-scaled diffusion coefficients
                Cxpmat=-Cxpmat(2:end-1,:);
                Cypmat=-Cypmat(:,2:end-1);

%                 barbeta=sqrt(min([Cxpmat(:);Cypmat(:)])*max([Cxpmat(:);Cypmat(:)]));

                %%
                %% Construct La (A)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mx=M;
My=M;

offyseq=[-Cypmat,zeros(Mx,1)];
offyseq=offyseq(:);
offxseq=[-Cxpmat;zeros(1,My)];
offxseq=offxseq(:);
A=-spdiags([offyseq,offxseq,-Cdmat(:),[0;offxseq(1:end-1)],[zeros(Mx,1);offyseq(1:end-Mx)]],[-My,-1,0,1,My],Mx*My,Mx*My);

                [mg(Level).A]=A;
                [mg(Level).Cxpmat]=Cxpmat;
                [mg(Level).Cypmat]=Cypmat;
                [mg(Level).Cdmat]=Cdmat;
                %%
                e=ones(m,1);
                Me=h*spdiags([e 4*e e]/6,-1:1,m,m); % 1D mass matrix
                %interpolation operator in matrix form
                Pn=(1/2)*spdiags([e 2*e e],-2:0,m,m);
                mg(Level).P=kron(Pn(:,1:2:end-2),Pn(:,1:2:end-2));

                mg(Level).M=[]; E=speye(m^2); mg(Level).shift=0;
                switch(smtype)
                    case 'Jacobi-C'
                        %Thm 4.2 in https://www.cs.ubc.ca/~greif/Publications/hg2021.pdf
                        c=2/h^2;cm=1/h^2;Lam=4/h^2+lamj;
                        b1=1-2*(c-cm)/Lam;b2=1+2*c/Lam;
                        mg(Level).shift=(abs(b1)/b1+abs(b2)/b2)/(abs(b1)+abs(b2));
                        mg(Level).factor=abs(b1-b2)/(abs(b1)+abs(b2));
                    case 'Vanka'
                        %mg(Level).M=@(r) Vanka2D_e(r,h,lamj); % matrix free M*r
                        %matrix version for faster CPU times
                        z=h^2*lamj;%complex number
                        aa=(z^2 + 8*z + 14)/(z^3 + 12*z^2 + 44*z + 48);
                        bb=1/(z^2 + 8*z + 12);cc=2/(z^3 + 12*z^2 + 44*z + 48);
                        Me=spdiags(h*[sqrt(cc)*e (2*bb/sqrt(cc))*e sqrt(cc)*e]/2,-1:1,m,m);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %mg(Level).M=kron(Me,Me)+h^2*(aa-bb^2/cc)*E;%large matrix
                        mg(Level).M=Me; mg(Level).shift=h^2*(aa-bb^2/cc);%small matrix
                    case 'Mass-shift'
                        z=h^2*lamj;
                        aa=(z + 2)/(z^2 + 4*z + 3);bb=1/(z^2 + 4*z + 3);
                        Me=h*spdiags([bb*e 2*aa*e bb*e]/2,-1:1,m,m); %shifted 1D mass matrix
                        mg(Level).M=kron(Me,Me);
                    case 'Vanka-e'
                        mg(Level).M=(3*kron(Me,Me)+h^2*E)/8; %Me in (23)
                    case 'Vanka-e2'
                        % Me2=h*spdiags([e 3*e e],-1:1,m,m); %
                        %mg(Level).M=(kron(Me2,Me2)+h^2*5*E)/28;
                        Me3=h*spdiags([e (10/3)*e e],-1:1,m,m); %optimal
                        mg(Level).M=(3*kron(Me3,Me3)+h^2*(32/3)*E)/80;
                    case 'Mass'
                        mg(Level).M=kron(Me,Me); %Mfe in (22) s
                end
            end
            A=mg(fineL).A;
     
            %%
             N=2^fineL;
             hx=(xmax-xmin)/N;
            x=rand(size(b));
             h=hx;
             r0=h*norm(b-A*x);rk=r0;
             err0=[]; it=1; err0(it)=rk; maxit=100;
            %x=full_mg_2d(mg,x,b,fineL,1,1,smtype);%Full Multigrid, faster
            %%
            
            %x=rand(size(b));
            %h=hx;
            %r0=h*norm(b-A*x);rk=r0;
           % err0=[]; it=1; err0(it)=rk;maxit=100;
            %x=full_mg_2d(mg,x,b,fineL,1,1,smtype);%Full Multigrid, faster
            %while((it<=maxit)&&(rk>=r0*tol))
             %   it=it+1;
              %  x=mg_iter_2d(mg,x,b,fineL,mu1,mu2,smtype);
               % rk=h*norm(b-A*x);
                %err0(it)=rk;
            %end
           % iter=length(err0)-1;%total V-cycle iterations
            %use AA-MG
            %[x,iter,err0]=AndAcc(@(x)mg_iter_2d(mg,x,b,fineL,mu1,mu2,smtype),x,10,maxit,tol,tol);
            %res=err0(end)/r0;
            
            %%
            
            while(it<=maxit)
                it=it+1;
                x=mg_iter_2d(mg,x,b,fineL,mu1,mu2,smtype);
                %rk=h*norm(b-A*x);            
            end
           
        end
       
    end
end

