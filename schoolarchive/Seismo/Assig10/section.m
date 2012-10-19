function section(seis,beg,dt,aflag);

% SECTION(SEIS,BEG,DT,AFLAG) plots a seismogram section of 
% seismograms in 2D array SEIS. BEG is begin time of section, DT
% is sample interval and AFLAG determines amplitude scaling. If
% AFLAG < 0 each trace is scaled to its maximum amplitude, if
% AFLAG = 0 each trace is scaled to the maximum amplitude of the
% entire section, and if AFLAG > 0 then each trace is scaled to
% that amplitude (note this option allows direct comparison of 
% plots produced by different calls to SECTION).
 
ny=size(seis,1);
nt=size(seis,2);
if nargin < 4
  aflag=-1;
end

ic=fix(rand(1)*7)+1;
col=['y','m','c','r','g','b','k'];

yaxe=[1:ny];
if aflag < 0
  for iy=1:ny
    xymat(iy,:)=iy-0.7*seis(iy,:)/max(abs(seis(iy,:))+0.0000001);
  end
elseif aflag == 0
  for iy=1:ny
    xymat(iy,:)=iy-8.0*seis(iy,:)/max(max(abs(seis))+0.0000001);
  end
else
  for iy=1:ny
    xymat(iy,:)=iy-1.0*seis(iy,:)/aflag;
  end
end
time=[0:nt-1]*dt+beg;
%plot(time,xymat,[col(ic),'-']);
plot(time,xymat,['k-']);
axis('ij');
xlab=xlabel('Time [s]');
ylab=ylabel('Trace #');
set(gca,'FontName','Helvetica','FontSize',14,'Clipping','off','layer','top');
set(xlab,'FontName','Helvetica','FontSize',14);
set(ylab,'FontName','Helvetica','FontSize',14);

