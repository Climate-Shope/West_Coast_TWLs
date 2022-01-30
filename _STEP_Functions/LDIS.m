function L=LDIS(T,h)

% L = LDIS(T,h)
% dispersion relationship
% Li Erikson USGS
%


g=9.81;
[Nr,Nc]=size(T);
[Nrh,Nch]=size(T);
if Nr ~= Nrh || Nc~=Nch
    disp('error - T and h need to be the same size array or matrix')
    return
end

for jj=1:Nr
    for kk=1:Nc
        omega=2*pi./T(jj,kk);
        D=omega.^2*h(jj,kk)/g;
        iter=0;
        iterm=50;
        error=1;

        if D>1
            xo=D;
        else
            xo=sqrt(D);
        end

        while(error)>0.001 && iter<10;
            F=xo-D/tanh(xo);
            DF=1+D/sinh(xo)^2;
            x1=xo-F/DF;
            error=abs((x1-xo)/xo);
            xo=x1;
            iter=iter+1;
        end

        if iter>iterm
            '10 iterations have been exceeded'
        else
            L(jj,kk)=2*pi*h(jj,kk)/x1;
        end
    end
end
