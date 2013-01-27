%Lab8
close all
clear all

run quakes;
save A;
vel=6;

E1=zeros(100,100);
for x=1:size(E1,1)
    for y=1:size(E1,2)
        res=0;
        T=mean(tq1-(((x-xcrd).^2+(y-ycrd).^2).^(1/2))/vel);
        for k=1:length(xcrd)
        res=res+(tq1(k)-(sqrt((x-xcrd(k))^2+(y-ycrd(k))^2)/vel)-T)^2;
        end
        E1(x,y)=res;
    end
end

E2=zeros(100,100);
for x=1:size(E2,1)
    for y=1:size(E2,2)
        res=0;
        T=mean(tq2-(((x-xcrd).^2+(y-ycrd).^2).^(1/2))/vel);
        for k=1:length(xcrd)
        res=res+(tq2(k)-(sqrt((x-xcrd(k))^2+(y-ycrd(k))^2)/vel)-T)^2;
        end
        E2(x,y)=res;
    end
end




[X(1),Y(1)]=find(E1==min(min(E1)))
[X(2),Y(2)]=find(E2==min(min(E2)))

%%Estimate uncertainty
%Compute sigma square
i=0;
T=mean(tq1-(((X(1)-xcrd).^2+(X(1)-ycrd).^2).^(1/2))/vel);
  for k=1:length(xcrd)
  i=i+(tq1(k)-(sqrt((X(1)-xcrd(k))^2+(Y(1)-ycrd(k))^2)/vel)-T)^2;
  end
sigma1=i/10;
i=0;
T=mean(tq2-(((X(2)-xcrd).^2+(Y(2)-ycrd).^2).^(1/2))/vel);
  for k=1:length(xcrd)
  i=i+(tq2(k)-(sqrt((X(2)-xcrd(k))^2+(Y(2)-ycrd(k))^2)/vel)-T)^2;
  end
sigma2=i/10;

%Kri square

K1=zeros(100,100);
for x=1:size(K1,1)
    for y=1:size(K1,2)
        res=0;
        T=mean(tq1-(((x-xcrd).^2+(y-ycrd).^2).^(1/2))/vel);
        for k=1:length(xcrd)
        res=res+(tq1(k)-(sqrt((x-xcrd(k))^2+(y-ycrd(k))^2)/vel)-T)^2;
        end
        K1(x,y)=res/sigma1;
    end
end


K2=zeros(100,100);
for x=1:size(K2,1)
    for y=1:size(K2,2)
        res=0;
        T=mean(tq2-(((x-xcrd).^2+(y-ycrd).^2).^(1/2))/vel);
        for k=1:length(xcrd)
        res=res+(tq2(k)-(sqrt((x-xcrd(k))^2+(y-ycrd(k))^2)/vel)-T)^2/sigma2;
        end
        K2(x,y)=res;
    end
end
figure (1)
imagesc(E1);
hold on
contour (K1,[18.31 18.31],'k')
plot(xcrd,ycrd,'o')


figure (2)
imagesc(E2);
hold on
contour (K2,[18.31 18.31],'k')
plot(xcrd,ycrd,'o')