function i_scal = convert_back(i_vec, deg_f)
% function i_scal = convert_back(i_vec, deg_f)

i_scal =0;

for j = 1:length(deg_f) 
    i_scal = i_scal + (i_vec(j)) * prod( deg_f(1:(j-1)));
end

end