function y=pad(x,n)

%returns column vector with the values of vector x padded (or truncated) to
% n elements

y=x(:);
sx=size(y);
sx=sx(1);

if sx < n
	y=[y.',zeros(1,n-sx)].';
end

y=y(1:n);

end 
