function s=powers_FFT(x,pow,flag)
% function s=powers_FFT(x,pow,flag)
% x = complex[1:n,1:2*m+1] OR complex{1:n}[1:2*m+1],
% n=number of vectors to multiply, 
% m=number of Fourier coefficients to consider,
% pow = integer[n]
% flag   a Bolean, default=1. If flag=1 the result s has the same size as
% x, otherwise it has maximum size. 
% this function calculates the product
% \Conv_{i=1}^{size(x)} x_i^{power(i)}
% were the power is considered to be a convolution.
% this result is obtained by transposing x in the real space, make a
% product there and going back to the Fourier space.

global use_intlab

if nargin<3
    flag=1;
end

if iscell(x)
    x_non_cell=zeros(max(size(x)),max(size(x{1})));
    for i=1:max(size(x))
        x_non_cell(i,:)=x{i};
    end
    x=x_non_cell;
end


%transpose=0;

if size(x,1) ~= size(pow,1)
    if size(x,2)==size(pow,1)
        x=x.';
    elseif size(x,1)==size(pow,2)
        pow=pow.';
    elseif size(x,2)==size(pow,2)
        pow=pow.';
        x=x.';
    else
        error('dimensions must agree')
    end
end

% x(i,:) has size 2m+1
m=length(x(1,:));
m=(m-1)/2;

if any(pow<0)
    s=0*x(1,:);
    s(m+1)=1;
    s=check_transpose(s,x);
    return
    %error('Powers cannot be negative')
elseif all(pow==0)
    s=0*x(1,:);
    s(m+1)=1;
    s=check_transpose(s,x);
    return
elseif sum(pow)==1
    s=x(pow==1,:);
    s=check_transpose(s,x);
    return
end



%if m~=floor(m)
%error('the size of x is not odd, there is a problem')
%end

% ta=zeros(size(x,2)+2*sum(pow)*m,1);
% if use_intlab
%     ta=intval(ta)';
% end
% tu=ta;

prod=1;
if use_intlab
    prod=intval(1);
end


% x(i,:) has size 2m+1
m=length(x(1,:));
m=(m-1)/2;
delta = ( floor( sum( pow ) /2) +1) *m;
len = length( x( 1,: ) );
cols = size( x,1 );
next_pow = ceil( log2( len +2*delta));
diff_len =  2 ^ next_pow - len ;

length_added_is_even=0;

if mod(diff_len,2)==1
    
    m=(diff_len-1)/2;
    y=[zeros(m+1,cols);x.';zeros(m,cols)].';
else
    
    length_added_is_even=1;

    m=(diff_len)/2;
    y=[zeros(m,cols);x.';zeros(m,cols)].';
end

for i=1:size(x,1)
    if pow(i)==0
        continue
    end
    ta=y(i,:);
%     if use_intlab
%         %ta=[intval(zeros(delta,1));x(i,:).';intval(zeros(delta,1))].';
%     else
%         ta=[zeros(delta,1);x(i,:).';zeros(delta,1)].';
%     end
    tu=(verifyfft_in(ta,-1));%ifftshift(ta(i,:)),-1));
    if use_intlab
        prod=prod.*(tu.^intval(pow(i)));
    else
        %     tu=ifft(ifftshift(ta));
        if 1==0%i==1
            prod=tu.^pow(1);
        else
            prod=prod.*(tu.^pow(i));
        end
    end
end
%if use_intlab
F=verifyfft_in(prod);
%else
%    F=fftshift(fft(prod));
%end


l=length(F);


if flag==1
    %if use_intlab
    %    s=((intval(4)*intval(m+1)+intval(1))^(sum(intval(pow))-intval(1))*F(delta+1:end-delta));
    %else
    if length_added_is_even
        s=(l^(sum(pow)-1)*F(m+1:end-m,:));
    else
        s=(l^(sum(pow)-1)*F(m+2:end-m,:));
    end
    %end
else
    
    %if use_intlab
    %    s=((intval(4)*intval(m+1)+intval(1))^(sum(intval(pow))-intval(1))*F(m+1:end-m));
    %else
        s=(l^(sum((pow))-(1))*F(2:end));
    %end
    
end

s=check_transpose(s,x);

return
