clear
clc

% make sure L is larger enough for converge
L=128;
kmin=1;
kmax=5;
theta=(1:L)./L*2*pi;
xx=cos(theta);
yy=sin(theta);
b=0.2;
c=0.2;
delta=pi/3;

%%%%%%%%%%%%%%共形变换z=z
z=xx*1+yy*sqrt(-1);
g=(1+2*b*z+3*c*exp(sqrt(-1)*delta)*z.^2)./sqrt(1+2*b^2+3*c^2);%求导
w=(z+b*z.^2+c*exp(sqrt(-1)*delta)*z.^3)./sqrt(1+2*b^2+3*c^2);%映射
uu=real(w);
vv=imag(w);
for i=1:L
    for j=1:L
        rho(i,j)=abs(w(i)-w(j));
    end
end

for i=1:L
    for j=1:L
        M(i,j)=(conj(z(i))*z(j)*conj(g(i))*g(j))/(abs(g(i)*g(j)))-1;
        N(i,j)=sqrt(-1)*((conj(z(i))*conj(g(i))*(w(i)-w(j)))/(abs(g(i))*abs(w(i)-w(j)))+(z(j)*g(j)*conj(w(i)-w(j)))/(abs(g(j))*abs(w(i)-w(j))));
        gg(i,j)=sqrt(abs(g(i)*g(j)));
    end
end

k=kmin:0.01:kmax;
for i=1:length(k)
    T=-pi*k(i)/(2*L).*gg.*(sqrt(-1).*M.*besselh(0,1,k(i).*rho)+N.*besselh(1,1,k(i)*rho));
    for i1=1:L
        T(i1,i1)=1;
    end
    D(i)=det(T);
    disp(i/length(k))
end


DD(1,:)=abs(D);
DD(2,:)=k;
load('africaDEk_2_0_500.mat')
plot(k,DD(1,:));hold on;plot(Ek(1:7),0,'r*')