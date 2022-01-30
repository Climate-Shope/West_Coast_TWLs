function c=hessigev(func,x,data),
ep=0.0001;
eps=ep*x;
n=length(x);
m=zeros(n,n);
for i=1:n,
    for j=1:n,
        x1=x;
        x1(i)=x1(i)+eps(i);
        x1(j)=x1(j)+eps(j);
        x2=x;
        x2(i)=x2(i)+eps(i);
        x2(j)=x2(j)-eps(j);
        x3=x;
        x3(i)=x3(i)-eps(i);
        x3(j)=x3(j)+eps(j);
        x4=x;
        x4(i)=x4(i)-eps(i);
        x4(j)=x4(j)-eps(j);
        
        
        m(i,j)=eval(['(' func '(x1,data)-' func '(x2,data)-' func '(x3,data)+' func '(x4,data)' ')' '/(4*eps(' num2str(i) ')*eps(' num2str(j) '))' ]);
    end
end

        c=inv(m);
        
		
