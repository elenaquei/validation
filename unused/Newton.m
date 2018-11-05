function [mu,T,x,y] = Newton( F, DF, mu, T, x, y)
tol = 10^-12;
max_iter = 10;
i=1;
RES(i) = norm(F(mu,T,x,y));
while i <= max_iter% && norm(F(mu,T,x,y))>tol
    
    vec_old=[mu,T,horiz(x),horiz(y)].';
    
    DF_findiff=fin_dif(F,mu,T,x,y);
    
    %problem = DF_findiff - DF(mu, T, x,y);
    
    vec_new = vec_old + DF_findiff\F(mu,T,x,y);
    
    mu = vec_new(1);
    T=vec_new(2);
    xy = vec_new(3:end);
    x=xy(1:length(xy)/2);
    y=xy(length(xy)/2+1:end);
    
    i=i+1;
    RES(i) = norm(F(mu,T,x,y));
end

res = norm(F(mu,T,x,y));
fprintf('res = %e\n',res)

return
end

function DF_findif = fin_dif(G,mu,T,x,y)
fun=@(x) G(x(1),x(2),x(3:2+length(y)),x(3+length(y):end));
diff = 10^-12;
vec = [mu,T,horiz(x),horiz(y)].';
e=eye(length(vec));

DF_findif=zeros(length(vec));

for k=1:length(vec)
        DF_findif(:,k) =( fun(vec+e(:,k)*diff) - fun(vec-e(:,k)*diff))/(2*diff);
end

end