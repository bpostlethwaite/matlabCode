function [Xp,Tp] = getdxdt(depth,alpha,p)


Xp = zeros(length(p),1);
Tp = zeros(length(p),1);


for jj = 1:length(p)
    
    ii = 1; % Set depth and velocity back to 0 for next p
    irtr = NaN;
    
    while ii <= length(depth)-1 && irtr ~= 2
        [dx,dt,irtr] = layerxt(p(jj),depth(ii+1) - depth(ii),...
            alpha(ii),alpha(ii+1));
        Xp(jj) = Xp(jj) + dx;
        Tp(jj) = Tp(jj) + dt;
        ii = ii+1;
    end
end

end