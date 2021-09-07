function [boundary_X,boundary_Y,Total_points]=Boundary(coeffi,L)
    a=coeffi(1);
    b=coeffi(2);
    delta=coeffi(3);
    Norm=sqrt(1+2*a^2+3*b^2);% normalization
    d_theta=2*pi/L;%rhe angle of two points
    theta=0:d_theta:(2*pi-d_theta);
    Z=exp(1i*theta);
    Z_polar=(Z+a*Z.^2+b*Z.^3*exp(sqrt(-1)*delta))/Norm;
    boundary_X=1*real(Z_polar);
    boundary_Y=1*imag(Z_polar);
    Total_points=length(boundary_X);
end