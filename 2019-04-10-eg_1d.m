clear all; 
close all; 
clc;


n = 32;
mres = 16;          % Number of sub-sampling total = n*mres


 
%% Creates Uniform Grid for Reconstruction
k = -n/2:n/2; k=k';
j=0:n-1;
x=-1+j/(n/2);
fx=zeros(n,1); 



for j=1:n
    if x(j)<0
        fx(j)=-x(j)-1;
    else if x(j)>=0
            fx(j)=-x(j)+1;
        end
    end
end



f_hat=zeros(n+1,1);

%Define Partition
j=0:mres*n-1;
quadx=-1+j/(mres*n/2);
f_quadx=zeros(mres*n,1);   
h=quadx(2)-quadx(1);    %delta x
for j=1:mres*n
    if quadx(j)<0
        f_quadx(j)=-quadx(j)-1;
    else if quadx(j)>=0
            f_quadx(j)=-quadx(j)+1;
        end
    end
end

for i=1:n+1
    f_hat(i)=0;
    for j=1:mres*n
        f_hat(i)=f_hat(i)+f_quadx(j)*exp(-1i*k(i)*pi*quadx(j))*h;
    end
end
f_hat=.5*f_hat;

S_n=zeros(n,n+1);
for i=1:n
    for j=1:n+1
        S_n(i,j)=exp(1i*k(j)*pi*x(i));
    end
end


alpha=32;
p=4;
G=zeros(n,n+1);
for i=1:n+1
    for j=1:n
        G(j,i)= exp(-alpha*(abs(k(i)/(n/2))).^p);
    end
end
F=G.*S_n;
f_filter=F*f_hat;
f_rec=S_n*f_hat;


h = figure; 
plot(x,fx,'ob',quadx(find(quadx<0)),f_quadx(find(quadx<0)),'b',quadx(find(quadx>=0)),f_quadx(find(quadx>=0)),'b','LineWidth',1,'MarkerSize',5)
set(gca,'FontSize',16)
axis([-1 1 -1.2 1.2])
h_legend = legend('Test Function','Uniform Grid');
set(h_legend,'FontSize',12,'Location','NorthWest')


h = figure; 
plot(x,real(f_rec),'ob',quadx(find(quadx<0)),f_quadx(find(quadx<0)),'b',quadx(find(quadx>=0)),f_quadx(find(quadx>=0)),'b','LineWidth',1,'MarkerSize',5)
set(gca,'FontSize',16)
axis([-1 1 -1.2 1.2])
h_legend = legend('Fourier Reconstruction','Test Function');
set(h_legend,'FontSize',12,'Location','NorthWest')


h = figure; 
plot(x,real(f_filter),'ob',quadx(find(quadx<0)),f_quadx(find(quadx<0)),'b',quadx(find(quadx>=0)),f_quadx(find(quadx>=0)),'b','LineWidth',1,'MarkerSize',5)
set(gca,'FontSize',16)
axis([-1 1 -1.2 1.2])
h_legend = legend('Filtered Fourier Reconstruction','Test Function');
set(h_legend,'FontSize',12,'Location','NorthWest')




