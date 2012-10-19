function [prof_lat prof_lon] = getprofiles(lat0,lon0,mainellipse,ring,offset,az1,az2,planet,units,num_pts,rng)

%% NOW THAT I HAVE REVERSED IT I NEED TO SWITCH THE ARC AROUND. ie, It
% is going from crater ellipse to fourth ring, thus fourth ring needs more
% points and arcs while major ring becomes the inside line.

for ii = 1:length(az2(:,1))
[elat,elon] = ellipse1(lat0,lon0,mainellipse,offset,az1(ii,:),planet,units,num_pts(1)); %an arc of Large ellipse lat and lon 
[elat2,elon2] = ellipse1(lat0,lon0,ring,offset,az2(ii,:),planet,units,num_pts(2)); % an arc Small ellipse lat and lon

Pt1 = [elat(round(median(length(az1(:,1))))).*ones(length(elat2),1),elon((round(median(length(az1(:,1)))))).*ones(length(elon2),1)]; %Pt1 is the middle lat and lon coord of the arc computed
% for the minor ellipse, copied over the same number of columns as the
% major ellipse arc
Pt2 = [elat2,elon2]; %Pt2 is the arc computed for the ring above

[dist,AZ] = distance(Pt1,Pt2,planet); % This computes the distance between the two arcs as well as the azimuthal angle
% between points on each arc 


[prof_lat(:,ii),prof_lon(:,ii)] = track1(elat(dist==min(dist)),elon(dist==min(dist)),AZ(dist==min(dist)),rng); %This finds the 
% lon and lat coord on the mainellipse arc which minimizes the distance
% between the two arcs, ie the norm, it uses the angle between this min
% distance arc to arc lon lat coords to draw (continue) a great circle line
% It stores the lat and lon coords.

end

end