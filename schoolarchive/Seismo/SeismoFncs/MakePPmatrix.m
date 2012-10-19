function [rtpp,rtppl,A]=MakePPmatrix(vp1,vs1,den1,vp2,vs2,den2,pp)
%%% Compute the reflection coefficient of the downgoing P wave
%%% You will obtain both the linearized results and the exact results
%%% Remember that A itself is linearized.
l=length(pp);
rtpp=zeros(l,16);
rtppl=zeros(l,8);
A=zeros(l,3);
for i=1:l
    rtpp(i,:)=RTCOEF(vp1,vs1,den1,vp2,vs2,den2,pp(i));
    [rtppl(i,:),A(i,:)]=RTCOEF_L(vp1,vs1,den1,vp2,vs2,den2,pp(i));
end
rtpp=rtpp(:,1);
rtppl=rtppl(:,1);