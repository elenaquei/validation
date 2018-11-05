function [mu,rho] = Newton_rho( F, DF, mu, rho)
tol = 10^-12;
max_iter = 10;
i=1;
RES(i) = norm(F(mu,rho));
while i <= max_iter% && norm(F(mu,T,x,y))>tol
    
    vec_old=[mu,horiz(rho)].';
    
    DF_findiff=fin_dif(F,mu,rho);
    
    problem = DF_findiff - DF(mu, rho);
    DF_mat=DF(mu,rho);
    
    vec_new = vec_old - DF_mat\F(mu,rho);
    
    mu = vec_new(1);
    rho=vec_new(2:end);
    
    i=i+1;
    RES(i) = norm(F(mu,rho));
end

res = norm(F(mu,rho));
fprintf('res = %e\n',res)

return
end

function DF_findif = fin_dif(G,mu,rho)
fun=@(x) G(x(1),x(2:end));
diff = 10^-12;
vec = [mu,horiz(rho)].';
e=eye(length(vec));

DF_findif=zeros(length(vec));

for k=1:length(vec)
        DF_findif(:,k) =( fun(vec+e(:,k)*diff) - fun(vec-e(:,k)*diff))/(2*diff);
end

end