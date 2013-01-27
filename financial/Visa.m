%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visa.m, Fun little plotting utility for the Visa balance
%
% Ben Postlethwaite       01/08/2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = 0;
date = '00/00/0000';
files = dir();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Visa Total and Dates out of file cell array
for ii = 1:length(files)
    d = files(ii).name;
    if length(d) > 3
        if strcmp(d(end-3:end),'.csv')
            fprintf('Processing %s\n',d)
            info = importdata(d);
            data = [data;info.data(:,3)];
            date = [date;cell2mat(info.textdata(:,1))];
            
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First values are zeros since I used array appending method
date(1,:) = [];
data(1) = [];
% Turn dates into datenumbers
xdata = datenum(date);

% Mean the heck out of the data
mdata = data;
for ii = 1:20
mdata = runmean(mdata,10);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT

figure1=figure(1);
plot(xdata,data,'r:',xdata,mdata)
datetick('x')
legend('Real Visa Balance','Running Mean of Balance','Location','NorthWest')
xlabel('Date')
ylabel('Balance in Dollars')
title(sprintf('The Visa Rollercoaster\nPart 1) The Climb'))
% Create textarrow
annotation(figure1,'textarrow',[0.402255639097744 0.485902255639098],...
    [0.640909814323607 0.461538461538462],'TextEdgeColor','none',...
    'String',{'Lunches, dinners and',' shopping spree.'});

% Create textarrow
annotation(figure1,'textarrow',[0.703007518796991 0.68609022556391],...
    [0.351785145888594 0.50748502994012],'TextEdgeColor','none',...
    'String',{'Ben Discovers','Granville Island'});

% Create textarrow
annotation(figure1,'textarrow',[0.803571428571428 0.773496240601504],...
    [0.587859416445623 0.708222811671087],'TextEdgeColor','none',...
    'String',{'Till the New Years,','fuck it.'});

% Create textarrow
annotation(figure1,'textarrow',[0.581766917293233 0.602443609022556],...
    [0.754968169761273 0.583554376657825],'TextEdgeColor','none',...
    'String',{'Cards go into','book-of-blood'});

