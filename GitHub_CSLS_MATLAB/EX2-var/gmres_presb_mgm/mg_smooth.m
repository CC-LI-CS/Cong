function [x]=mg_smooth(A,x0,b,nv,smoother,M,shift)
x=x0;
switch smoother
    case {'Mass','Mass-shift'} %mass smoother
        w=3/4;
        for k = 1:nv
            x=x + w*M*(b - A*x);
        end
    case 'Vanka'
        %w=24/25;
        r=shift; w=(20*r^2+20*r+48)/(2*r^2+21*r+50);
        for k = 1:nv
            res=(b - A*x);
            %x=x + w*M*(b - A*x); %avoid construct large M matrix
            Me_res=M*reshape(res,sqrt(length(b)),[])*M;
            x=x+w*(Me_res(:)+shift*res);
        end
    case 'Vanka-e'
        %w=24/25;
        r=shift; w=(20*r^2+20*r+48)/(2*r^2+21*r+50);
        for k = 1:nv
            x=x + w*M*(b - A*x);
        end
    case 'Vanka-e2'
        w=0.5253; %w=0.5801;
        for k = 1:nv
            x=x + w*M*(b - A*x);
        end
    case 'Jacobi'
        w=4/5;
        dinv = 1./diag(A);
        for k = 1:nv
            x=x + w*dinv.*(b - A*x);
        end
    case 'Jacobi-C'
        w=shift;
        dinv = 1./diag(A);
        for k = 1:nv
            x=x + w*dinv.*(b - A*x);
        end
    case 'GS' %G-S Smoother
        L = tril(A);
        for k = 1:nv
            x=x + (L\(b - A*x));
        end
    case 'RB-GS'
        dinv = 1./diag(A);
        for k = 1:nv
            x(1:2:end)=x(1:2:end) + dinv(1:2:end).*(b(1:2:end) - A(1:2:end,:)*x);
            x(2:2:end)=x(2:2:end) + dinv(2:2:end).*(b(2:2:end) - A(2:2:end,:)*x);
        end
end
end