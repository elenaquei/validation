classdef coefs
    properties %(SetAccess=private)
        size_scalar
        size_vector
        deg_scalar
        deg_vector
        non_zero_el
        powers_scalars
        powers_vectors
        value
        %linear_coefs %2-cell array for linear part % NOT YET assigned in constructors
    end
    methods
        % COEFS
        function alpha=coefs(Sscal,Svec,Dscal,Dvec,Non0,Pscal,Pvec,val)
            % set the number of scalars and vectors, the degrees of each,
            % the number of non zero coefficients, the powers of the scalars
            % and vectors and the values connected to them
            % Pscal, Pvec and val must be cell array of dimension Svec
            if nargin >=2
                if Sscal<0
                    error('number of scalars must be non-negative')
                elseif Svec<=0
                    error('number of vectors must be positive')
                end
                alpha.size_scalar=Sscal;
                alpha.size_vector=Svec;
            end
            if nargin >= 4
                alpha.deg_scalar=Dscal;
                alpha.deg_vector=Dvec;
            end
            if nargin ==8
                
                if max(size(Non0)) < Svec
                    error('Non-zero coeffs must be at least the number of vectors')
                elseif size(Pscal{1},1)~= Sscal
                    if size(Pscal{1},2)==Sscal
                        Pscal=Pscal.';
                    else
                        error('dimensions do not agree')
                    end
                elseif size(Pvec) ~= Svec
                    error('dimensions do not agree')
                    %elseif max(size(val)) ~= Svec
                    %    error('dimensions do not agree')
                end
                alpha.non_zero_el=Non0;
                alpha.powers_scalars=Pscal;
                alpha.powers_vectors=Pvec;
                alpha.value=val;
            end
        end
        % COEFS
        
        % COPY_COEFS
        function beta=copy_coefs(alpha)
        % beta is in all aspects fo copy of alpha
        beta.size_scalar=alpha.size_scalar;
        beta.size_vector=alpha.size_vector;
        beta.deg_scalar=alpha.deg_scalar;
        beta.deg_vector=alpha.deg_vector;
        beta.non_zero_el=alpha.non_zero_el;
        beta.powers_scalars=alpha.powers_scalars;
        beta.powers_vectors=alpha.powers_vectors;
        beta.value=alpha.value;
        
        end
        % COPY_COEFS
        
        % EQ
        function bol= eq(alpha, beta)
           bol=0;
           if alpha.size_scalar ~= beta.size_scalar
               return
           elseif alpha.size_vector ~= beta.size_vector
               return
           else
               for i=1:alpha.size_vector
                   if any(alpha.powers_scalars{i}~=beta.powers_scalars{i})
                       return
                   elseif any(alpha.powers_vectors{i}~= beta.powers_vectors{i})
                       return
                   elseif any(alpha.value{i}~= beta.value{i})
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
                alpha.size_scalar;
                alpha.size_vector;
                alpha.deg_scalar;
                alpha.deg_vector;
                alpha.non_zero_el;
                alpha.powers_scalars;
                alpha.powers_vectors;
                alpha.value;
            catch
                return
            end
            bool=1;
            return
        end
        % ISCOEF
    end
end
