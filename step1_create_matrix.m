function [G_all, D_all, in]=step1_create_matrix(InSAR_overlap, meanll, gps_h, gps_3d)
% Step 1 - generating a jumbo G martix  

% INPUTS:
% data ------- data matrix of three overlapping InSAR tracks [n*15] [lon1, lat1, v1, inc1, phi1, lon2, lat2,
% v2, inc2, phi2, lon3, lat3, v3, inc3, phi3] 
% meanll ----- central lon/lat [n*2]
% gps_h ------ data matrix of gps horizontal components and overlapping InSAR tracks [m*15]  [lon1, lat1, v1, inc1, phi1, lon2, lat2,
% v2, inc2, phi2, lon3, lat3, v3, inc3, phi3] 
% gps_3d ----- data matrix of gps 3d components and the overlapping InSAR
% track [p*10] [lon1, lat1, v1, inc1, phi1, lon2, lat2,
% v2, inc2, phi2] 

% OUTPUTS:
% G_all ------ jumbo G matrix
% D_all ------ all data points
% in --------- parameter index

% By Lin Shen -- University of Leeds


%%%Senario 1 - InSAR points in the overlapping region that having both ascending and descending measurements

v_match_121_filter=InSAR_overlap;

v_match_121_filter(:,1:2)=v_match_121_filter(:,1:2)-repmat(meanll,length(v_match_121_filter),1);
v_match_121_filter(:,6:7)=v_match_121_filter(:,6:7)-repmat(meanll,length(v_match_121_filter),1);
v_match_121_filter(:,11:12)=v_match_121_filter(:,11:12)-repmat(meanll,length(v_match_121_filter),1);

Da11=v_match_121_filter(:,3);
Ga11=[-1*cosd(v_match_121_filter(:,5)).*sind(v_match_121_filter(:,4)),sqrt(sind(v_match_121_filter(:,5)).^2.*(sind(v_match_121_filter(:,4)).^2)+cosd(v_match_121_filter(:,4)).^2),ones(length(Da11),1),v_match_121_filter(:,1),v_match_121_filter(:,2),zeros(length(Da11),6)];

Da21=v_match_121_filter(:,8);
Ga21=[-1*cosd(v_match_121_filter(:,10)).*sind(v_match_121_filter(:,9)),sqrt(sind(v_match_121_filter(:,10)).^2.*(sind(v_match_121_filter(:,9))).^2+cosd(v_match_121_filter(:,9)).^2),zeros(length(Da11),3),ones(length(Da11),1),v_match_121_filter(:,6),v_match_121_filter(:,7),zeros(length(Da11),3)];

Dd121=v_match_121_filter(:,13);
Gd121=[cosd(v_match_121_filter(:,15)).*sind(v_match_121_filter(:,14)),sqrt(sind(v_match_121_filter(:,15)).^2.*(sind(v_match_121_filter(:,14))).^2+cosd(v_match_121_filter(:,14)).^2),zeros(length(Da11),6),ones(length(Da11),1),v_match_121_filter(:,11),v_match_121_filter(:,12)];


G_all1=[Ga11;Ga21;Gd121];

G_add_gps=G_all1(:,3:end);
G_no_gps=G_all1(:,1:2);

G_p1=zeros(length(G_no_gps)/3,length(G_no_gps)/3*2);
for i=1:size(G_p1,1)
    for j=i:i
    G_p1(i,(j-1)*2+1:j*2)=G_no_gps(i,:);
    end
end

G_p2=zeros(length(G_no_gps)/3,length(G_no_gps)/3*2);
for i=1:size(G_p2,1)
    for j=i:i
    G_p2(i,(j-1)*2+1:j*2)=G_no_gps(i+size(G_p2,1),:);
    end
end

G_p3=zeros(length(G_no_gps)/3,length(G_no_gps)/3*2);
for i=1:size(G_p3,1)
    for j=i:i
    G_p3(i,(j-1)*2+1:j*2)=G_no_gps(i+size(G_p3,1)*2,:);
    end
end

G_temp=[G_p1;G_p2;G_p3];

D_temp=[Da11;Da21;Dd121];


G_insar=[G_temp,G_add_gps];
D_insar=[D_temp];


%%%%%%%Senario 2 - InSAR points that having ascending and descending vaules, and overlapped with GPS data that has the horizontal component only 

gps_h(:,1:2)=gps_h(:,1:2)-repmat(meanll,size(gps_h,1),1);
gps_h(:,6:7)=gps_h(:,6:7)-repmat(meanll,size(gps_h,1),1);
gps_h(:,11:12)=gps_h(:,11:12)-repmat(meanll,size(gps_h,1),1);

Ggps_h_e1=[ones(size(gps_h,1),1),zeros(size(gps_h,1),1),zeros(size(gps_h,1),9)];
Gsar_asc1=[-cosd(gps_h(:,10)).*sind(gps_h(:,9)),sqrt(sind(gps_h(:,10)).^2.*(sind(gps_h(:,9))).^2+cosd(gps_h(:,9)).^2),zeros(size(gps_h,1),3),ones(size(gps_h,1),1),gps_h(:,6),gps_h(:,7),zeros(size(gps_h,1),3)];
Gsar_desc1=[cosd(gps_h(:,15)).*sind(gps_h(:,14)),sqrt(sind(gps_h(:,15)).^2.*(sind(gps_h(:,14))).^2+cosd(gps_h(:,14)).^2),zeros(size(gps_h,1),6),ones(size(gps_h,1),1),gps_h(:,11),gps_h(:,12)];


Dgps_h_e1=gps_h(:,3);
Dsar_asc1=gps_h(:,8);
Dsar_desc1=gps_h(:,13);


% 
Ggps_h_all=[Ggps_h_e1;Gsar_asc1;Gsar_desc1];
Dgps_h_all=[Dgps_h_e1;Dsar_asc1;Dsar_desc1];

G_gps_1=Ggps_h_all(:,1:2);
G_gps_2=Ggps_h_all(:,3:end);
%  

G_p1=zeros(length(G_gps_1)/3,length(G_gps_1)/3*2);
for i=1:size(G_p1,1)
    for j=i:i
    G_p1(i,(j-1)*2+1:j*2)=G_gps_1(i,:);
    end
end

G_p2=zeros(length(G_gps_1)/3,length(G_gps_1)/3*2);
for i=1:size(G_p2,1)
    for j=i:i
    G_p2(i,(j-1)*2+1:j*2)=G_gps_1(i+size(G_p2,1),:);
    end
end
G_p3=zeros(length(G_gps_1)/3,length(G_gps_1)/3*2);
for i=1:size(G_p3,1)
    for j=i:i
    G_p3(i,(j-1)*2+1:j*2)=G_gps_1(i+size(G_p3,1)*2,:);
    end
end

G_gps_n=[G_p1;G_p2;G_p3];

G_gps_1=G_gps_n;
G_gps_h=[G_gps_2,G_gps_1];
D_gps_h=Dgps_h_all;


%%%%%%%Senario 3 - InSAR points that having ascending or descending vaule, and overlapped with GPS data that has the horizontal component plus the vertical 


gps_3d(:,1:2)=gps_3d(:,1:2)-repmat(meanll,size(gps_3d,1),1);
gps_3d(:,6:7)=gps_3d(:,6:7)-repmat(meanll,size(gps_3d,1),1);

Ggps_3d_e1=[ones(size(gps_3d,1),1),zeros(size(gps_3d,1),1),zeros(size(gps_3d,1),9)];
Ggps_3d_u1=[zeros(size(gps_3d,1),1),sqrt(sind(gps_3d(:,10)).^2.*(sind(gps_3d(:,9))).^2+cosd(gps_3d(:,9)).^2),zeros(size(gps_3d,1),9)];
Gsar_asc1=[-cosd(gps_3d(:,10)).*sind(gps_3d(:,9)),sqrt(sind(gps_3d(:,10)).^2.*(sind(gps_3d(:,9))).^2+cosd(gps_3d(:,9)).^2),ones(size(gps_3d,1),1),gps_3d(:,6),gps_3d(:,7),zeros(size(gps_3d,1),6)];
Dgps_3d_e1=gps_3d(:,3);
Dgps_3d_u1=gps_3d(:,4).*sind(gps_3d(:,10)).*sind(gps_3d(:,9))+gps_3d(:,5).*cosd(gps_3d(:,9));
Dsar_asc1=gps_3d(:,8);

Ggps_3d_all=[Ggps_3d_e1;Ggps_3d_u1;Gsar_asc1];
Dgps_3d_all=[Dgps_3d_e1;Dgps_3d_u1;Dsar_asc1];

G_gps_1=Ggps_3d_all(:,1:2);
G_gps_2=Ggps_3d_all(:,3:end);
 
G_p1=zeros(length(G_gps_1)/3,length(G_gps_1)/3*2);
for i=1:size(G_p1,1)
    for j=i:i
    G_p1(i,(j-1)*2+1:j*2)=G_gps_1(i,:);
    end
end

G_p2=zeros(length(G_gps_1)/3,length(G_gps_1)/3*2);
for i=1:size(G_p2,1)
    for j=i:i
    G_p2(i,(j-1)*2+1:j*2)=G_gps_1(i+size(G_p2,1),:);
    end
end
G_p3=zeros(length(G_gps_1)/3,length(G_gps_1)/3*2);
for i=1:size(G_p3,1)
    for j=i:i
    G_p3(i,(j-1)*2+1:j*2)=G_gps_1(i+size(G_p3,1)*2,:);
    end
end

G_gps_n=[G_p1;G_p2;G_p3];

G_gps_1=G_gps_n;
G_gps_3d=[G_gps_2,G_gps_1];
D_gps_3d=Dgps_3d_all;

%%%%%%%%%%%Combine matrices of three senarios to a jumbo one

ss=size(G_gps_3d,1);
ww=size(G_gps_3d,2);

G_gps_all_tmp=G_gps_h;
G_gps_all_tmp(size(G_gps_h,1)+1:size(G_gps_h,1)+ss,size(G_gps_h,2)+1:size(G_gps_h,2)+ww-9)=G_gps_3d(:,9+1:end);
G_gps_all_tmp(size(G_gps_h,1)+1:size(G_gps_h,1)+ss,1:9)=G_gps_3d(:,1:9);
G_gps_all=G_gps_all_tmp;

D_gps_all=[D_gps_h;D_gps_3d];


G_all=G_insar;
G_all(size(G_insar,1)+1:size(G_insar,1)+size(G_gps_all,1),size(G_insar,2)-9+1:size(G_insar,2)+size(G_gps_all,2)-9)=G_gps_all;
D_all=[D_insar;D_gps_all];
in=size(G_insar,2)-9+1;

