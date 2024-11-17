function y = ilu_vec(Cxpmat,Cypmat,Cdmat,eigB,tv,M)

Mx=M;
My=M;

offyseq=[-Cypmat,zeros(Mx,1)];
offyseq=offyseq(:);
offxseq=[-Cxpmat;zeros(1,My)];
offxseq=offxseq(:);
K=-spdiags([offyseq,offxseq,-Cdmat(:),[0;offxseq(1:end-1)],[zeros(Mx,1);offyseq(1:end-Mx)]],[-My,-1,0,1,My],Mx*My,Mx*My);
I=speye(size(K));
setup.type='nofill';  
[L,U]=ilu( K+eigB*I ,setup);
y=U\(L\tv);
yy=(K+eigB*I)\tv;
ny=norm(y-yy);
end