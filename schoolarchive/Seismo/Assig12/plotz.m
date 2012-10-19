close all

for ii = 1:length(p)
    hold on
    %plot(w,W{ii}(1))
    a(ii) = W{ii}(1);
end

plot(w,a)