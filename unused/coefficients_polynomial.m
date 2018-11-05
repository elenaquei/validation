classdef coefficients_polynomial
    properties %(SetAccess=private)
        size_vector_field % M 
        deg_f   % vector of length M
        max_deg  % integer
        non_zero_el % vector of lenght M
        powers  % cell {M}(non_zero_el,M)
        powers_lambda  % cell {M}(non_zero_el)
        coefficients % cell {M}(non_zero_el)
        phase_condition % vector {M}          - phase condition: y(0) * p0 = 0
        amplitude_condition % vector {M}      - amplitude condition: y(0) * p1 = 1
    end
    methods
        % COEFS
        function alpha=coefficients_polynomial(dim, deg, non_zero, power, power_lambda, coef, p0, p1)
            
                if length(dim) ~=1
                    error('dimension of the vector field must be scalar');
                end
                alpha.size_vector_field = dim;
            if nargin>1
                if length(deg)==1
                    deg= deg*ones(dim,1);
                end
                alpha.deg_f=deg;
                if nargin>2
                    if numel(non_zero)~=dim
                        error('number of elements must a vector of length dim');
                    end
                    alpha.non_zero_el = non_zero;
                    if nargin>3
                        if ~iscell(power)
                            error('power must be a cell of dimension dim , dim x non_zero_el')
                        end
                        for i =1:dim
                            if any(size(power{i})~= [non_zero(i),dim])
                                error('Dimensions do not match')
                            end
                        end
                        alpha.powers = power;
                        if nargin>4
                            if ~iscell(power_lambda)
                                error('power must be a cell of dimension dim , dim x non_zero_el')
                            end
                            for i =1:dim
                                if any(length(power_lambda{i})~= non_zero)
                                    error('Dimensions do not match')
                                end
                            end
                            alpha.powers_lambda = power_lambda;
                            if nargin>5
                                if ~iscell(coef)
                                    error('coef must be a cell of dimension dim, non_zero_el')
                                end
                                for i =1:dim
                                    if length(coef{i})~= non_zero
                                        error('Dimensions do not match')
                                    end
                                end
                                alpha.coefficients = coef;
                                if nargin >6
                                    if min(size(p0))~= 1
                                        error('p0 must be a vector')
                                    elseif max(size(p0)) ~= dim
                                        error('length of p0 bust be equal to dimesion of the vector field')
                                    end
                                    alpha.phase_condition = p0;
                                    if nargin>7
                                    if min(size(p1))~= 1
                                        error('p1 must be a vector')
                                    elseif max(size(p1)) ~= dim
                                        error('length of p1 bust be equal to dimesion of the vector field')
                                    end
                                    alpha.amplitude_condition = p1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
        end
        % COEFS
        
        % SCALAR_CONDITION
        function beta = scalar_condition(alpha, p0, p1)
        %function beta = scalar_condition(alpha, p0, p1)
        % adds the scalar condititons to alpha
        
        beta = alpha;
        beta.amplitude_condition = p1;
        beta.phase_condition = p0;
        end
        % SCALAR_CONDITION
        
        
        
        % COPY_COEFS
        function beta=copy_coefs(alpha)
        % beta is in all aspects a copy of alpha
        
        beta.size_vector_field =alpha.size_vector_field; 
        beta.deg_f = alpha.deg_f  ;
        beta.max_deg = alpha.max_deg; 
        beta.non_zero_el = alpha.non_zero_el ;
        beta.powers = alpha.powers;
        beta.powers_lambda = alpha.powers_lambda;
        beta.coefficients = alpha.coefficients;
        beta.phase_condition = alpha.phase_condition ;
        beta.amplitude_condition = alpha.amplitude_condition ;
        
        end
        % COPY_COEFS
        
        % EQ
        function bol= eq(alpha, beta)
           bol=0;
           if beta.size_vector_field ~=alpha.size_vector_field 
               return
           elseif any(beta.deg_f ~= alpha.deg_f )
               return
           elseif beta.max_deg ~= alpha.max_deg
               return
           elseif any(alpha.phase_condition ~= beta.phase_condition)
               return
           elseif any(alpha.amplitude_condition ~= beta.amplitude_condition)
               return
           else
               for i=1:alpha.size_vector
                   if any(alpha.powers_lambda{i}~=beta.powers_lambda{i})
                       return
                   elseif any(alpha.powers{i}~= beta.powers{i})
                       return
                   elseif any(alpha.coefficients{i}~= beta.coefficients{i})
                       return
                   elseif alpha.non_zero_el{i} ~= beta.non_zero_el{i}
                       return
                   end
               end
           end
           bol=1;
        end
        % EQ
        
        % ISCOEF
        function bool = iscoef(alpha)
            bool=0;
            try
                alpha.size_vector_field;
                alpha.deg_f;
                alpha.max_deg;
                alpha.non_zero_el;
                alpha.powers;
                alpha.powers_lambda;
                alpha.coefficients;
                alpha.phase_condition;
                alpha.amplitude_condition;
            catch
                return
            end
            bool=1;
            return
        end
        % ISCOEF
    end
end
