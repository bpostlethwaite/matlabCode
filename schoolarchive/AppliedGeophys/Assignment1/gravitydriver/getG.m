function [G] = getG(x,y,z,x0,y0,dx,dy,dz)

h = dx*dy*dz;

g=zeros(length(x),length(y),length(z));
G=zeros(length(x)*length(y)*length(z),length(x0)*length(y0));
count = 1;
for kk = 1:length(y0)
for mm = 1:length(x0)   




        for jj = 1:length(z)
                for ii = 1:length(y)
                   for qq = 1:length(x)
                        g(qq,ii,jj) =(z(jj)*h)*1/(((x(qq)-x0(mm))^2 + (y(ii) - y0(kk))^2 + z(jj)^2).^(3/2));
            
                   end
                end
        end
G(:,count) = reshape(g,[length(x)*length(y)*length(z),1]);
count = count +1;       
        
end
end

end