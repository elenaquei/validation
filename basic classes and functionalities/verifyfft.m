function Z = verifyfft(z,sign)
%VERIFYFFT    Verified 1-dimensional FFT
%
%   res = verifyfft(z,sign)
%
%   z     input vector or matrix
%         length of z must be a power of 2 
%   sign   1 forward FFT (default)
%         -1 inverse FFT
% 
%As in Matlab, the inverse FFT is scaled such that forward and inverse FFT
%are inverse operations.
%For matrix input, FFT is performed on each column; row vector input
%is converted into column vector. 
%For N-dimensional FFT apply verifyfft N times.
% 

% written  09/24/14     S.M. Rump  (based no Marcio Gameiro's code)
%
  
  [n,col] = size(z);
  if n==1
    if col==1
      Z = intval(z);
      return
    else
      z = z(:);
      n=col;
      col=1;
    end
  end
    
  if nargin==1
    sign = 1;       % default: forward
  end
  
  % check dimension
  log2n = round(log2(n));
  if 2^log2n~=n
    error('length must be power of 2')
  end
  
  % bit-reversal
  % v = bin2dec(fliplr(dec2bin(0:n-1,log2n))) + 1
  f = 2^(log2n-1);
  v = [0;f]; 
  for k=1:log2n-1
    f = 0.5*f;
    v = [ v ; f+v ];
  end
  z = z(v+1,:);
  
  % Danielson-Lanczos algorithm
  Z = intval(z);
  Index = reshape(1:n*col,n,col);
  nmax = 16384;     % maximum in fft_data
  if n<0%<=nmax
    load fft_data   % roots of unity in  r +/- d
    Phi = midrad(r(1:nmax/n:nmax),d);
    if sign==-1
      Phi = (Phi.')';      
    end
  else
    % compute roots of unity, division exact because n is power of 2
    theta = intval('pi') * ( sign*(0:(n-1))'/n ); 
    Phi = cos(theta) + 1i*sin(theta);
  end
  v = 1:2:n;
  w = 2:2:n;
  t = Z(w,:);
  Z(w,:) = Z(v,:) - t;
  Z(v,:) = Z(v,:) + t;
  
  for index=1:(log2n-1)     % Executed log2(n) times
    m = 2^index;
    m2 = 2*m;
    vw = reshape(1:n,m2,n/m2);
    v = vw(1:m,:);
    w = vw(m+1:m2,:);
%     t = bsxfun(@times,exp(1i*pi*(0:m-1)'/m),Z(w));  % doesn't work for intervals
%     theta = intval('pi') * (sign*(0:(m-1))'/m);     % division exact because m=2^p
%     t = exp(1i*theta) .* Z(w);
    indexv = reshape(Index(v(:),:),m,col*n/m2);
    indexw = reshape(Index(w(:),:),m,col*n/m2);
%     t = repmat(Phi(1:n/m:end),1,n/m2*col);
    t = Phi(1:n/m:end,ones(1,n/m2*col)) .* Z(indexw);   % Tony's trick
    Z(indexw) = Z(indexv) - t;
    Z(indexv) = Z(indexv) + t;
  end
  
  Z = [Z(1,:); flipud(Z(2:end,:))];
  if sign==-1
    Z = Z/n;        % error-free since n is a power of 2
  end
  
end
