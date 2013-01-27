function out = getequal(f1,f2,w,modes)

out = NaN(1,modes);
f = f1 - real(f2);
count = 1;
for ii = 1:length(f) -1
    if f(ii) < 0 && f(ii+1) > 0
        out(count) = (w(ii) + w(ii+1)) / 2;
        count = count + 1;
    end   
end
 
end