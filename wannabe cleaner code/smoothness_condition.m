function flag = smoothness_condition(x0, x_dot_0, x1,x_dot_1, r)
global nu

scal_norm_1 = 1 + 0*x0.scalar;

abs_K = abs(-x0.nodes:x0.nodes);

vec_norm_1 =nu .^-abs_K;

Xi_vec_norm_1 = Xi_vector(scal_norm_1, repmat(vec_norm_1,x0.size_vector,1));

vec_in_ball = Xi_vec2vec(interval_Xi(x0,x1) + interval_Xi(-r* Xi_vec_norm_1,r*Xi_vec_norm_1));

vec_x0 = Xi_vec2vec(x0);
vec_x1 = Xi_vec2vec(x1);

x_dot_Delta = x_dot_1-x_dot_0;

neq0 = vec_in_ball(:).'*x_dot_Delta(:)+vec_x0(:).'*x_dot_0(:)-vec_x1(:).'*x_dot_1(:);

if real(neq0)>0 || real(neq0)<0 || imag(neq0)>0 || imag(neq0)<0 
    flag = 1;
else
    flag = 0;
end

return