clear all
close all

n = 1000;
s = 0:1/(n-1):1;
x = s;
p = 1:10:1000;

d = 1;

for ii = 1:length(s)
    Kst(ii,:) = d./((d.^2 + (s(ii) - x).^2).^3/2);
end

for ii = 1:length(p)
    fpt(ii,:) = sin(2*pi*p(ii)*x);
end


for jj = 1:length(p)
    for ii = 1:length(s)
        gps(jj,ii) = sum(Kst(ii,:).*fpt(jj,:));
        I(jj,ii) = sum(x.*fpt(jj,:));
    end
end



gph = char('b','r','g','k');
for ii = 1:length(p)
    hold on
    figure(1)
    plot(x,fpt(ii,:),sprintf('%s',gph(ii)))
    figure(2)
    plot(s,gps(ii,:),sprintf('%s',gph(ii)))
end
figure(3)
plot(p,I)