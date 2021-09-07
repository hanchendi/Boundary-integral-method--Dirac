function  [M,N,length_2,S,M_A,N_A,Ealpha_S,kai1,kai2,alpha_S,alpha,NX,NY,TX,TY,Kai_1,Kai_2]=Structure(boundary_X,boundary_Y,Total_points)
format long
TN5=Total_points;
X=boundary_X;
Y=boundary_Y;
XX=zeros(1,TN5+2);
YY=zeros(1,TN5+2);
XX(1)=X(TN5);
YY(1)=Y(TN5);
XX(2+TN5)=X(1);
YY(2+TN5)=Y(1);
XX(2:TN5+1)=X(1:TN5);
YY(2:TN5+1)=Y(1:TN5);
for g=2:TN5+1 
    SQL_1=(XX(g+1)-XX(g))^2+(YY(g+1)-YY(g))^2;%前项差分
    SQL_2=(XX(g)-XX(g-1))^2+(YY(g)-YY(g-1))^2;%后项差分
    SQL=(XX(g+1)-XX(g-1))^2+(YY(g+1)-YY(g-1))^2;%前后
    normL=sqrt(SQL);
    normL_1=sqrt(SQL_1);
    normL_2=sqrt(SQL_2);
    S(g-1)=(normL_1+normL_2)/2;
    TX(g-1)=(XX(g+1)-XX(g-1))/normL;TY(g-1)=(YY(g+1)-YY(g-1))/normL;%切线 of Sｐｏｉｎｔ
    NX(g-1)=(YY(g+1)-YY(g-1))/normL;NY(g-1)=-(XX(g+1)-XX(g-1))/normL;%normal of S point
    %NX(g-1)=-(YY(g+1)-YY(g-1))/normL;NY(g-1)=(XX(g+1)-XX(g-1))/normL;%normal of S point
end

K=TN5 ;                      %  relation between i and j
ones_K=ones(K,1);
X_2=ones_K*X;
Y_2=ones_K*Y;
X_X=-X_2+X_2';%%X_X(m,n)=x(m)-x(n);
Y_Y=-Y_2+Y_2';%%Y_Y(m,n)=y(m)-y(n)
X_Y=X_X+1i*Y_Y;
clear X_2 Y_2
length_2=abs(X_Y);%矩阵长的行或列
SX=X_X./length_2;%%SX(m,n)=x(m)-x(n);
SY=Y_Y./length_2;
clear X_X Y_Y

for g=1:K 
    SX(g,g)=TX(g);              %two point direction   
    SY(g,g)=TY(g);
end
%     SNX=-SY;
%     SNY=SX;

alpha=zeros(1,K);
for g=1:K                     %calculate alpha
     slope=NY(g)/ NX(g);
  if  NX(g)>=0&& NY(g)>=0
      alpha(g)=atan(slope);
  end
  if  NX(g)<0&& NY(g)>=0
      alpha(g)=atan(slope)+pi;
  end
  if  NX(g)<0&& NY(g)<0
  %    alpha(g)=atan(slope)-pi;
  alpha(g)=atan(slope)+pi;
  end
  if   NX(g)>=0&&NY(g)<0
    %   alpha(g)=atan(slope);
   alpha(g)=atan(slope)+2*pi;
  end
end
% define kai
kai1=zeros(K,K);kai2=zeros(K,K);

 NX_S=ones_K*NX;NX_S=NX_S';%%NX_S(m,n)=NS(m);
 NY_S=ones_K*NY;NY_S=NY_S';
 TX_S=ones_K*TX;TX_S=TX_S';%%TX_S(m,n)=TX(m);
 TY_S=ones_K*TY;TY_S=TY_S';
 J=NX_S.*SX+NY_S.*SY;%%J(m,n)=NX_S(m)*(SX(m)-SX(n))+NY_S(m)*(SY(m)-SY(n))
 J=J';
 JNSN=TX_S.*SX+TY_S.*SY;%%JNSN(m,n)=TX_S(m)*(SX(m)-SX(n))+TY_S(m)*(SY(m)-SY(n));
 JNSN=JNSN';
% clear NX_S NY_S
%tic
 alpha_S=ones_K*alpha;
 %alpha_S=-alpha_S+alpha_S';%more right
 alpha_S=-alpha_S+alpha_S';%alpha_S(m,n)=alpha(m)-alpha(n);
 alpha_S=alpha_S';
sign_JNSN=sign(JNSN);
one_A=ones(K,K);
kai1=sign_JNSN.*acos(J)+pi*(one_A-sign_JNSN);%(0,2*pi)
kai2=kai1+alpha_S;%kai plus kai2(m,n)=kai(m,n)+alpha(m)-alpha(n);alpha_S(m,n)=alpha(m)-alpha(n);
 
   
    S_S=ones_K*S;%S_S a line include all S (s1,s2 ,s3,...,sn)
   % S_S=S_S';%%S_S a column include all S (s1,s2 ,s3,...,sn)
     %alpha_S=alpha_S';%%%change
    % alpha_S=-alpha_S;
   Ealpha_S=exp(-1i*alpha_S);
   M_A=(1/4)*(Ealpha_S-1);
   M=(1/4)*(Ealpha_S-1).*S_S;
   Kai_1=exp(1i*kai1);
   Kai_2=exp(-1i*kai2);
   N_A=(1i*1/4)*(Kai_1+Kai_2);
   N_A=1/2*(N_A+N_A');
   N=N_A.*S_S;
  % N=(1i*1/4)*(Kai_1+Kai_2).*S_S;
  % clear Kai_1 Kai_2
   
    for n=1:K
         M(n,n)=0;
         N(n,n)=0;
    end
    
end  