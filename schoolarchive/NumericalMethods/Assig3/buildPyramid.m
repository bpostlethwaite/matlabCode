function [Z] = buildPyramid(N,h)

for ii = 1:N
    for jj = 1:N
        Z(ii,jj) = h * abs(  (min([ii,jj,abs(N-ii),abs(N-jj)]) - 1)   );
    end
end

end