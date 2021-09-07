clear
clc

% make sure L is larger enough for converge

load('Dirac_african_B.mat');
mode=2;
k=Dirac_african_B(mode);
b=0.2;
c=0.2;
delta=pi/3;
L=300;
r_out=0.98;
Norm=sqrt(1+2*b^2+3*c^2);% normalization
d_theta=2*pi/L;%rhe angle of two points
theta=0:d_theta:(2*pi-d_theta);
Z=exp(1i*theta);
Z_polar=(Z+b*Z.^2+c*Z.^3*exp(sqrt(-1)*delta))/Norm;
boundary_X=1*real(Z_polar);
boundary_Y=1*imag(Z_polar);

[phi_1_L,phi_2_L,coordinate_X,coordinate_Y,normalized]=eig_state_di(boundary_X,boundary_Y,k,b,c,delta,L,r_out);
  
current=conj(phi_1_L).*phi_2_L;
u1=real(current);
u2=imag(current);
figure(1)
mesh(coordinate_X,coordinate_Y,abs(phi_1_L).^2+abs(phi_2_L).^2);
figure(2)
T=4;
plot(coordinate_X(L,:),coordinate_Y(L,:),'--k');hold on
quiver(coordinate_X(L,1:T:L),coordinate_Y(L,1:T:L),u1(L,1:T:L),u2(L,1:T:L),'b','linewidth',1.5);
axis([-1 1.5 -1.2 1])
axis off

%save([pwd,'/Compare/psi1_BX_',num2str(mode),'.mat'], 'phi_1_L');
%save([pwd,'/Compare/psi2_BX_',num2str(mode),'.mat'], 'phi_2_L');
save([pwd,'/psi1_Boundary_',num2str(mode),'.mat'], 'phi_1_L');
save([pwd,'/psi2_Boundary_',num2str(mode),'.mat'], 'phi_2_L');