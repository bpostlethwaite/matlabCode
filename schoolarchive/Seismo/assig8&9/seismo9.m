% Assignment 9
% Ben Postlethwaite
clear all; close all
warning off

colour = {'vw','+k','*r','xk'};
for j=1:2
    a = [0,2];
    a = a(j);
    seismo8;
    %send variable to second implementation of prior code so I don't
    %      overwrite figures
    param  = [1,6]; % The stability mechanism for reducing step size
    steps = [12,50]; % Two step sizes for the two different regimes
    for p = 1:2
        mGuess = [20, 20  
                  20, 80 
                  10, 10];
        for m = 1:2
            for ii = 1:steps(j)              
                figure(p+a)
                hold on
                plot(mGuess(1,m),mGuess(2,m),sprintf('%s',colour{m}))

                r = tq{p} - sqrt((mGuess(1,m)-xcrd).^2+(mGuess(2,m)-ycrd).^2)/vel - mGuess(3,m);
                g = sqrt((xcrd - mGuess(1,m)).^2 + (ycrd-mGuess(2,m)).^2)*vel/param(j);
                G = [(mGuess(1,m)-xcrd)./(g) ,(mGuess(2,m)-ycrd)./(g) ,ones(length(ycrd),1)];
                
                mcor = G\r;
                mGuess(:,m) = mGuess(:,m) + mcor;
            end
        end
    end
end
