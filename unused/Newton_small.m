function [x] = Newton_small( F, DF,x)
tol = 10^-12;
max_iter = 10;
i=1;
RES(i) = norm(F(x));
while i <= max_iter% && norm(F(mu,T,x,y))>tol
    
    DF_findiff=fin_dif(F,x);
    
    problem = DF_findiff - DF(x);
    DF_mat=DF(x);
    
    x = x - DF_mat\F(x);
    
    i=i+1;
    RES(i) = norm(F(x));
end

res = norm(F(x));
fprintf('res = %e\n',res)

return
end


function DF_findif = fin_dif(G,x)
fun=@(x) G(x);
diff = 10^-12;
e=eye(length(x));

DF_findif=zeros(length(x));

for k=1:length(x)
        DF_findif(:,k) =( fun(vert(x)+e(:,k)*diff) - fun(vert(x)-e(:,k)*diff))/(2*diff);
end
end