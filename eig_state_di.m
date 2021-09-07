function [phi_1_L,phi_2_L,coordinate_X,coordinate_Y,normalized]=eig_state_di(boundary_X,boundary_Y,k,b,c,delta,L,r_out)
format long

minus=0.5*10^(-2);
nx=L;
ny=L;
rr=linspace(0,r_out,L);
rr=rr';
theta=linspace(0,2*pi,L);
dr=rr(2)-rr(1);
dtheta=theta(2)-theta(1);

xx1=rr*cos(theta);
yy1=rr*sin(theta);
z=xx1+sqrt(-1)*yy1;
w=(z+b*z.^2+c*exp(sqrt(-1)*delta)*z.^3)./sqrt(1+2*b^2+3*c^2);
g=(1+2*b*z+3*c*exp(sqrt(-1)*delta)*z.^2)./sqrt(1+2*b^2+3*c^2);
coordinate_X=real(w);
coordinate_Y=imag(w);
for i=1:L
    for j=1:L
        normalized(i,j)=rr(i)*abs(g(i,j))^2*dr*dtheta;
    end
end
Coor_line_X=coordinate_X(:);
Coor_line_Y=coordinate_Y(:);
%save([pwd,'/coordinate_X','_','.mat'], 'coordinate_X');
%save([pwd,'/coordinate_Y','_','.mat'], 'coordinate_Y');
%save([pwd,'/normalized.mat'], 'normalized');
save([pwd,'/g.mat'], 'g');
N_CXY=nx*ny;
N_Y=ny;
N_X=nx;
%[xv,yv,N_point_in_boundary]=Boundary(coeffi,2*boundary_points);
[M,N,length_2,S,M_A,N_A,Ealpha_S,kai,kai2,alpha_S,alpha,NX,NY,TX,TY,Kai_1,Kai_2]=Structure(boundary_X,boundary_Y,L);
   %disp(Ek)
   %disp(length_2)
   save([pwd,'/alpha','.mat'], 'alpha');
   [~,u_s]=matrix_D_us(k,L,length_2,M,N);

   boundary_points=L;
   S=S';
  normal_SS=u_s.*S;%[boundary,1]
  phi_1L=zeros(N_CXY,1,'double');
  phi_2L=zeros(N_CXY,1,'double');
  poss_0=zeros(N_CXY,1,'double');
  phi_1=zeros(ny,nx,'double');
  phi_2=zeros(ny,nx,'double');
  possibility=zeros(ny,nx,'double');
 alpha=alpha';
Normal_S_alpha=normal_SS.*exp(-1i*alpha);
ones_B=ones(1,boundary_points,'double');
ones_NY=ones(1,N_Y,'double');
boundary_XR=ones_NY'*boundary_X;
boundary_YR=ones_NY'*boundary_Y;
NXR=ones_NY'*NX;
NYR=ones_NY'*NY;
TXR=ones_NY'*TX;
TYR=ones_NY'*TY;

boundary_points=boundary_points;
ones_B=ones(1,boundary_points,'double');
ones_NXY=ones(N_CXY,1,'double');
boundary_XL=ones_NXY*boundary_X;
boundary_YL=ones_NXY*boundary_Y;
NLX=ones_NXY*NX;
NLY=ones_NXY*NY;
TLX=ones_NXY*TX;
TLY=ones_NXY*TY;


boundary_XL=ones_NXY*boundary_X;
boundary_YL=ones_NXY*boundary_Y;
%size_boundary_X=size(boundary_XL)

cood_xL=Coor_line_X*ones_B;
cood_yL=Coor_line_Y*ones_B;

R_X=boundary_XL-cood_xL;
R_Y=boundary_YL-cood_yL;
rou_L=R_X+1i*R_Y;
length_L=abs(rou_L);
rou_L_x=R_X./length_L;
rou_L_y=R_Y./length_L;
cos_kai_L=rou_L_x.*NLX+rou_L_y.*NLY;
sin_kai_L=rou_L_x.*TLX+rou_L_y.*TLY;
kai_L_1=cos_kai_L+1i*sin_kai_L;
kai_L_2=cos_kai_L-1i*sin_kai_L;
hankle_L_1=besselh(1,1,k*length_L);%besselh(0,1,length)
t=1;
for i=1:length(length_L(:,1))
    for j=1:length(length_L(1,:))
        if length_L(i,j)<= minus
            label(t,1)=i;
            label(t,2)=j;
            t=t+1;
        end
    end
end
save([pwd,'/label.mat'], 'label');
hankle_L_1(label(:,1),label(:,2))=0;
hankle_L_0=besselh(0,1,k*length_L);%besselh(1,1,length)
        hk1_L=hankle_L_1.*kai_L_1;%%[N_Y,boundary]
        hk2_L=hankle_L_1.*kai_L_2;
        phi_1L=(hk1_L-hankle_L_0)*normal_SS;
        phi_2L=(hk2_L+hankle_L_0)*Normal_S_alpha;
        phi_1L=(1i*k/4)*phi_1L;
        phi_2L=(k/4)*phi_2L;
  phi_1_L=reshape(phi_1L,ny,nx);   
  phi_2_L=reshape(phi_2L,ny,nx);
  phi_1_L=conj(phi_1_L);
  phi_2_L=conj(phi_2_L);
end