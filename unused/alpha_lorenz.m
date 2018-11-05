function alpha_new=alpha_lorenz(alpha_input,pho)
% function alpha_new=alpha_lorenz(alpha,pho)
%
% INPUT
% alpha      coef for Lorenz system
% pho        new paramenter for Lorenz system
%
% OUTPUT
% alpha_new  coef for Lorenz system with pho as in input

if ~iscoef(alpha_input)
    error('first input must be the coefficients for Lorenz')
end
if any(size(pho)~=1)
    error('the second input must be a scalar')
end

alpha_new=alpha_input;
alpha_new.value{2}(2)=pho/(2*pi);
if alpha_input.size_scalar==2
    alpha_new.value{2}(2)=1/(2*pi);
end
alpha_coef=alpha_new;

if alpha_coef.size_scalar==1
    save('lorenz.mat','alpha_coef')
    save('saved elements/lorenz.mat','alpha_coef')
    save('lorentz.mat','alpha_coef')
    save('saved elements/lorentz.mat','alpha_coef')
elseif alpha_coef.size_scalar==2
    save('lorenz_cont.mat','alpha_coef')
    save('saved elements/lorenz_cont.mat','alpha_coef')
    save('lorentz_cont.mat','alpha_coef')
    save('saved elements/lorentz_cont.mat','alpha_coef')
else
    error('size_scalar too big')
end
return