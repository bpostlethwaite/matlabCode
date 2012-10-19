% this function takes a vector and zero_pad, or truncate, to length 2^13, or 8192

function Y = zero_pad(X)

if ( length(X) >= 2^13)
	Y = X(1:2^13);
else 
	Y = X;
	Y(length(X)+1:2^13) = zeros(1,2^13-length(X));
end
