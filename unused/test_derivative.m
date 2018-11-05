function DF= test_derivative(xBar)
% DF= test_derivative(xBar)
%

x=xBar.vector(1,:);
y=xBar.vector(2,:);
omega=xBar.scalar(1);
nodes=xBar.nodes;

oneN=diag(-nodes:nodes);
one = diag(1+0*[-nodes:nodes]);

if length(xBar.scalar)==1
    lambda=1;
    
    DF=[0                  ones(1,2*nodes+1)            zeros(1,2*nodes+1)
        y.'/(2*pi)         1i*oneN                      omega*lambda/(2*pi)*one
        -x.'/(2*pi)        -omega*lambda/(2*pi)*one     1i*oneN];
    
    if size(DF,1)~=size(DF,2)
        error('aarg')
    elseif size(DF,1)~= 1+(2*xBar.nodes+1)*xBar.size_vector
        error('2')
    end
    
else
    lambda=xBar.scalar(2);
    DF=[0                    0                     ones(1,2*nodes+1)            zeros(1,2*nodes+1)
        lambda*y.'/(2*pi)    omega*y.'/(2*pi)      1i*oneN                      omega*lambda/(2*pi)*one
        -lambda*x.'/(2*pi)   -omega*x.'/(2*pi)     -omega*lambda/(2*pi)*one     1i*oneN];
    
    if size(DF,1)+1~=size(DF,2)
        error('aarg')
    elseif size(DF,2)~= 2+(2*xBar.nodes+1)*xBar.size_vector
        error('2')
    end
end