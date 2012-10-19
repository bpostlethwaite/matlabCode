% Ben Postlethwaite
% Assignment 4 
% Data Analysis
%{
(1) The given data file (NLR.mat file in Matlab 7 format or text file)
contains two predictors x1 and x2 and predictand data y (80 observations).
Develop a nonlinear regression model of y as a function of x1 and x2.
Briefly describe the approach taken. Forecast ytest using the new test
predictor data x1test and x2test provided. The smaller is the RMSE between
your ytest and the true value, the higher will be your mark for this
exercise. [It would be most convenient if you also email me an Excel
    spreadsheet containing your 40 forecasted ytest values in a column
    to save me from having to type in your values to compute the RMSE].
%}
clear all
close all

load NLR.mat
data = [x1,x2];
targets = y';

load LowRMSindex.mat
addpath(genpath('/home/bpostlet/Dropbox/ComLinks/programming/matlab/DataAnalysis/Stochastic_Bosque'))

%% Data Division

% Divides data into training data, valuation and test data following an
% 60/20/20 rule
%[trainP,valP,testP,trainInd,valInd,testInd] = dividerand(data');
%[trainT,valT,testT] = divideind(t,trainInd,valInd,testInd);
iters = size(BetterINDtrain,1)*40;

h = waitbar(0,sprintf('completing %i iterations with three neural simulations, be patient',iters));
q = 0;
c = 1;
for ii = 1:iters
    waitbar(ii/iters,h)
    
    % 80% are more representative of the full data set
    % Randomize Index for valuation versus Training data
    %if iters <= 20
    %    intraining = LowRMSindex(1,:);
    %else
    %    intraining = LowRMSindex(2,:);
    %end
    %
    if mod(q,40) == 0
        M = length(targets);
        %m = round(0.8*M);
        %intraining = randperm(M);
        %intraining = sort(intraining(1:m));
        intraining = BetterINDtrain(c,:);
        notintraining = setdiff([1:M],intraining);
        
        %intrainINDEX(c,:) = intraining;
        %notintrainingINDEX(c,:) = notintraining;
        c = c+1;
    end
    q = q + 1;
    xindex(ii) = ceil(q/40); % make an index so we can select the appropriate
    % intrainINDEX of intraining values. ie xindex(342) = 54
    % intrainINDEX(54) = {values of index at desired low point}
    
    %}

    
    %% Bayseian Network
    % Map to [-1 1]
    [pn,ps] = mapminmax(data(intraining,:)');
    [tn,ts] = mapminmax(targets(intraining));
    
    bnet=newff(pn,tn,20,{},'trainbr');
    bnet.trainParam.showWindow = 0;
    bnet.divideFcn ='';
    bnet.performFcn = 'sse';
    bnet.trainParam.show = 10;
    bnet.trainParam.epochs = 50;
    randn('seed',192736547);
    bnet = init(bnet);
    [bnet,btr]=train(bnet,pn,tn);
    
    % Try out trained Baysian network on remaining data
    p2n = mapminmax('apply',data(notintraining,:)',ps);
    an = bnet(p2n);
    a = mapminmax('reverse',an,ts);
    
    error = targets(:,notintraining) - a;
    Bayrms(ii) = norm(error)/sqrt(numel(error));
    fprintf('Bayseian RMS error is %2.2f \n',Bayrms(ii))
    
    % Now run test data to get test y out
    p2n = mapminmax('apply',[x1test,x2test]',ps);
    an = bnet(p2n);
    ytestbays(ii,:) = mapminmax('reverse',an,ts);
    
 
    %% PLot
    %{
    figure(34)
    plot3(x1(notintraining),x2(notintraining),y(notintraining),'o',...
        x1(notintraining),x2(notintraining),a,'.');%,...
    %x1(notintraining),ones(length(x2(notintraining)),1),y(notintraining),...
    %ones(length(x1(notintraining)),1),x2(notintraining),y(notintraining))
    grid on
    axis square
    %}
    %% Feed Forward network
    %{
Early stopping. Set parameters to reduce convergence time.
For early stopping, you must be careful not to use an algorithm that
converges too rapidly. If you are using a fast algorithm (like trainlm),
set the training parameters so that the convergence is relatively slow.
For example, set mu to a relatively large value, such as 1,
and set mu_dec and mu_inc to values close to 1, such as 0.8 and 1.5,
respectively. The training functions trainscg and trainbr usually work
well with early stopping.
With early stopping, the choice of the validation set is also important.
The validation set should be representative of all points in the training set.
    %}
    %{
    net = feedforwardnet;
    net.trainFcn = 'trainscg';
    net.trainParam.showWindow = 0;
    net.divideFcn = 'divideInd';
    net.divideParam.trainInd = intraining;
    net.divideParam.valInd = notintraining(1:floor(0.5*length(notintraining)));
    net.divideParam.testInd = notintraining(ceil(0.5*length(notintraining)):end);
    net.trainParam.mu_dec = 0.8;
    net.trainParam.mu_inc = 1.5;
    [net,tr] = train(net,data',targets);
    
    a = bnet(data(notintraining,:)');
    error = targets(:,notintraining) - a;
    
    FeedForwardrms(ii) = norm(error)/sqrt(numel(error));
    fprintf('Feed Forward Conjugate Gradient RMS error is %2.2f \n',FeedForwardrms(ii))
    %}
    
    %{
yout = net(data);
trout = yout(tr.trainInd);
vout = yout(tr.valInd);
tsout = yout(tr.testInd);
trTarg = targets(tr.trainInd);
vTarg = targets(tr.valInd);
tsTarg = targets(tr.testInd);
plotregression(trTarg,trout,'Train',vTarg,vout,'Validation',...
tsTarg,tsout,'Testing')
    %}
    
    %% Random Forest
    %{
intraining = 1:60;
notintraining = 61:80;
    %}
    RandomForest = Stochastic_Bosque(data(intraining,:),targets(:,intraining),'ntrees',100);
    [output votes] = eval_Stochastic_Bosque(data(notintraining,:),RandomForest);
    ytestranf(:,ii) = eval_Stochastic_Bosque([x1test,x2test],RandomForest);
    error = targets(:,notintraining)'-output;
    RanForestrms(ii) = norm(error)/sqrt(numel(error));
    fprintf('RandomForest RMS error is %2.2f \n',RanForestrms(ii))
    %}
end
figure(534)
%plot([1:iters],Bayrms,[1:iters],FeedForwardrms,[1:iters],RanForestrms)
%legend('Baysian Neural Net','Feed Forward Neural Net','Random Forest')
%title('RMS error with Validation Data')
plot([1:iters],Bayrms,[1:iters],RanForestrms)
legend('Baysian Neural Net','Random Forest')
title('RMS error with Validation Data')

ytestbays = ytestbays';
save('NLR.mat','ytestbays','ytestranf','RanForestrms','Bayrms','-append')

figure(353)
%plot([1:iters],Bayrms + FeedForwardrms + RanForestrms)
plot([1:iters],Bayrms + RanForestrms)
%{
figure(2354)
plot3(x1(notintraining),x2(notintraining),y(notintraining),'o',...
    x1(notintraining),x2(notintraining),output,'.');%,...
    %x1(notintraining),ones(length(x2(notintraining)),1),y(notintraining),...
    %ones(length(x1(notintraining)),1),x2(notintraining),y(notintraining))
grid on
axis square
%}