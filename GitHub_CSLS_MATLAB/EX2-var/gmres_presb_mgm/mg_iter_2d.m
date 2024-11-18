function [x]=mg_iter_2d(mg,x0,b,level,pre,post,smoother)
%fprintf('+++++++++++++++++++Level[%d] Begin+++++++++++++++++++++\n',level);
A=mg(level).A;
global LEVEL0;
if(level==LEVEL0) %Coarest level
    x=A\b;
else
    % presmooth
    x = mg_smooth(A,x0,b,pre,smoother,mg(level).M,mg(level).shift);
    % Restrict residual
    r = b - A*x;
    %rc = mg_restricth(r,h);
    %rc = rest(r);
    rc=(mg(level).P)'*r/4;
    % coarse grid correction
    cc = mg_iter_2d(mg,zeros(size(rc)),rc,level-1,pre,post,smoother);
    cc = mg_iter_2d(mg,cc,rc,level-1,pre,post,smoother); %add this for W cycle
    %x = x + intp2(cc);
    %x=x + mg_interph(cc,h/2);
    x = x + mg(level).P*cc; %
    % postsmooth
    %fprintf('Level[%d],Before Postsmooth res=%1.2e\n',level,norm(b-A*x));
    x = mg_smooth(A,x,b,post,smoother,mg(level).M,mg(level).shift);
end
end
