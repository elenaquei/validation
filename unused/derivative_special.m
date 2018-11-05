function der = derivative_special( y1, y2, mu, a, omega, p0, p1)


one_long = ones(1,length(y1));
p01 = one_long*p0(1);
p02 = one_long*p0(2);
p11 = one_long*p1(1);
p12 = one_long*p1(2);

n = (length(y1) -1)/2;
K = -n:n;

y1 = vert(y1);
y2 = vert(y2);
y1_y1 = vert(conv(y1,y1,'same'));
y1_y2 = vert(conv(y1,y2,'same'));
y2_y2 = vert(conv(y2,y2,'same'));
yy1 = 0*y1;
yy1(n+1)=1;

 d5 =[   y1, -1i*vert(K).*y1, - 6*y1_y1*a-2*y2_y2*a, -4*y1_y2*a , toeplitz(-18*y1_y1*a^2-6*y2_y2*a^2)+diag(-1i*K*omega), toeplitz(-yy1-12*y1_y2*a^2)];
 d6= [  y2, -1i*vert(K).*y2, -4*y1_y2*a, -2*y1_y1*a-6*y2_y2*a, toeplitz(yy1-12*y1_y2*a^2), toeplitz(-6*y1_y1*a^2-18*y2_y2*a^2)+diag(-1i*K*omega)];
    
der = [ 0,0,0,0,p01,p02;
    0,0,0,0,p11,p12;
    0,0,0,-1,0*one_long,0*one_long;
    0,0,1,0,0*one_long, 0*one_long;
    -d5;
    -d6];



