clear all
ITER=300;
Name_system='vanderpol_cont';
for K=40:10:100
s=sprintf('./simulations/%s_validation_%d_%d.mat',Name_system,ITER,K);
load(s)
figure(1)
semilogy(K,max(delta_vec),'*');
title('max(delta)')
max(delta_vec)
hold on
figure(2)
plot(K,end_time,'*');
title('total time')
hold on
figure(3)
plot(K,iter,'*');
title('number of iterations')
hold on
figure(4)
semilogy(K,x1.scalar(2),'*');
title('final validated mu')
hold on
figure(5)
plot(delta_vec);
title('delta')
max(delta_vec)
hold on
end