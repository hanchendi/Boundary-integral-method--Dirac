function [D,u_s]=matrix_D_us(k,L,length,M,N)

   length_k=k*length;
   minus=0.5*10^(-2);
   t=1;
   for i=1:L
       for j=1:L
           if i~=j && length(i,j)<=minus
              label(t,1)=i;
              label(t,2)=j;
              t=t+1;
          end
      end
   end
   D=-k*(1i*M.*besselh(0,1,length_k)+N.*besselh(1,1,length_k));
   D(label(:,1),label(:,2))=-k*(1i*M(label(:,1),label(:,2)).*besselh(0,1,length_k(label(:,1),label(:,2))));
                      
   for n=1:L
       D(n,n)=1;
   end
   % D=conj(D');
   % det_k=det(D);det_k=abs(det_k);
   %[U,S,V]=svd(D);
   %u_s=V(:,size(D,2));
   [s11 s12]=eig(D);
   s13=diag(abs(s12));
   t=find(s13==min(s13));
   u_s=s11(:,t);  
   u_s0=u_s;
   %save([pwd,'/u_s0.mat'], 'u_s0');
end