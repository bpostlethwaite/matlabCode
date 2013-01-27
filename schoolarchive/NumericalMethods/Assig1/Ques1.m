clc
clear all
close all

h = 10.^-[20:-1:0];
x0 = 1.2;

f  = repmat(-sin(x0),length(h),1);
f1 = (cos(x0 + h) - cos(x0)) ./ h;
f2 = -(2*sin(x0 + h/2).*sin(h/2)) ./ h;

fprintf(repmat('%-16s',1,6),'h','real','f1','error in f1','f2','error in f2') %#ok<*CTPCT,*PRTCAL>
fprintf('\n')
for i = 1:length(h)
    fprintf([repmat('%+2.6e   ',1,6),' \n'],h(i),f(i),f1(i),abs(f(i)-f1(i)),f2(i),abs(f(i)-f2(i)))
end