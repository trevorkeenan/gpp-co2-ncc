% this script will...
% recreate Fig. 3
% from Keenan et al. 2023 https://doi.org/10.1038/s41558-023-01867-2


% mapping functions are provided by M_map:
% Pawlowicz, R., 2020. "M_Map: A mapping package for MATLAB", 
% version 1.4m, [Computer software], 
% available online at www.eoas.ubc.ca/~rich/map.html.


clear all
close all

saveFigures=0; % set to 1 to save figures

startLen = 10; % years
endLen = startLen-1;
mapMax = 125;

addpath('./functions')
addpath('./functions/m_map/')
printStatement='print(fig1,''-dpdf'', fname)';


% calculate the area of each latitude
lats=min(-89.5):0.5:max(90);
nlats = length(lats);
area=zeros(nlats,1);
earthellipsoid = almanac('earth','ellipsoid','m','sphere');
for ii=1:nlats
    lat1 = lats(ii);
    lat2 = lat1+1;
    area(ii) = areaquad(lat1,1,lat2,2,earthellipsoid); %m^2
end


%% 1. RS-GPP

filename=strcat('./data/dataFile.mat');
load(filename)
% conversions
gToPg=10^15;    % 10^15 gC = 1 PgC

lonlatland=lonlat;

% set the dimensions / lonlat of the entire grid
lat=-89.75:0.5:89.75;
lon=-179.75:0.5:179.75;
m=length(lat);
n=length(lon);

% calculate the area of each of the land pixels (here in m2)
lat1=(lonlatland(:,2)+0.25);
lat2=(lonlatland(:,2)-0.25);
lon1=((lonlatland(:,1)-0.25));
lon2=((lonlatland(:,1)+0.25));

earthellipsoid = referenceSphere('earth','m');
areaPR = areaquad(lat1,lon1,lat2,lon2,earthellipsoid); % area is in m2

for ii=1:size(i.GPP_MPI_CO2dyn,2)
    GPP_MPI_CO2dyn(:,ii)=i.GPP_MPI_CO2dyn(:,ii); % (gC m-2 yr-1) * (m-2)
    GPP_MOD17_CO2dyn(:,ii)=i.GPP_MOD17_CO2dyn(:,ii); % (gC m-2 yr-1) * (m-2)
    GPP_MPI_noCO2(:,ii)=i.GPP_MPI_noCO2(:,ii);
    GPP_MOD17_noCO2(:,ii)=i.GPP_MOD17_noCO2(:,ii);
end

for ii= 1:2
    switch ii
        case 1
            cGPP = GPP_MPI_CO2dyn;
        case 2
            cGPP = GPP_MOD17_CO2dyn;
    end
    startx=nanmean(cGPP(:,1:startLen),2);
    endx=nanmean(cGPP(:,end-endLen:end),2);
    diffx=(endx./startx);
    diffAbsRS(:,ii) = endx-startx;
    diffRelRS(:,ii) = (endx-startx)./startx;
    
    dataToMap=diffAbsRS(:,ii)';
    [ensembleDiffAbsRS_mapData, lon_b, lat_b] = regrid(lonlat,dataToMap);
    tmp = repmat(area,[1 720]);
    diffLat_RS(:,ii) = squeeze(nansum( tmp.*ensembleDiffAbsRS_mapData,2))./1e15;   % get latitudinal totals
    
end

diffLat_RSmean = nanmean(diffLat_RS,2);
diffLat_RSstd = nanstd(diffLat_RS,1,2);


clear startx endx
for ii= 1:2
    switch ii
        case 1
            cGPP = GPP_MPI_noCO2;
        case 2
            cGPP = GPP_MOD17_noCO2;
    end
    startx=nanmean(cGPP(:,1:startLen),2);
    endx=nanmean(cGPP(:,end-endLen:end),2);
    diffAbsRSnoC02(:,ii) = endx-startx;
    diffRelRSnoC02(:,ii) = (endx-startx)./startx;
    
    
    dataToMap=diffAbsRSnoC02(:,ii)';
    [ensembleDiffAbsRSnoCO2_mapData, lon_b, lat_b] = regrid(lonlat,dataToMap);
    tmp = repmat(area,[1 720]);
    diffLat_RSnoCO2(:,ii) = squeeze(nansum( tmp.*ensembleDiffAbsRSnoCO2_mapData,2))./1e15;   % get latitudinal totals
      
end

diffLat_RSnoCO2mean = nanmean(diffLat_RSnoCO2,2);
diffLat_RSnoCO2std = nanstd(diffLat_RSnoCO2,1,2);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot the latidudinal distribution


scrsz = get(0,'ScreenSize');
fig1 =figure('Position',[1 scrsz(4) scrsz(3)/3 scrsz(4)]);
hold on

colorx2={[0.8 0.7 0.6]*1.2,'k','r','r','r',[0.6,0.6,0.6]};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add the RS data (mean +- stdev)
scenario = 2;
y=flipud(diffLat_RSmean)*2;
E=flipud(diffLat_RSstd)*2; 

% truncate to get on  grid:-60:0.5:90
y=y(60:end);
E=E(60:end);

[l pH2]=boundedLine(y, 1:length(y), E,'orientation', 'horiz');
hnew = outlinebounds(l, pH2);
set(hnew,'color',colorx2{scenario})
set(l,'color',colorx2{scenario})
set(pH2,'FaceColor',colorx2{scenario})
set(pH2,'FaceAlpha',0.2)
set(pH2,'DisplayName','RS with CO_2')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add the RS data (mean +- stdev)
scenario = 3;
y=flipud(diffLat_RSnoCO2mean)*2;
E=flipud(diffLat_RSnoCO2std)*2; 

% truncate to get on  grid:-60:0.5:90
y=y(60:end);
E=E(60:end);

[l pH3]=boundedLine(y, 1:length(y), E,'orientation', 'horiz');
hnew = outlinebounds(l, pH3);
set(hnew,'color',colorx2{scenario})
set(l,'color',colorx2{scenario})
set(pH3,'FaceColor',colorx2{scenario})
set(pH3,'FaceAlpha',0.2)
set(pH3,'DisplayName','RS without CO_2')


% vertical lines
yLim=get(gca,'YLim');
x=0;
l1=line([x x],yLim);
set(l1,'LineStyle','--','Color','k')

xlim([-0.15 1.0])
ylim([0 300])

set(gca,'YTick',[10,60,110,160,210,260])
set(gca,'YTickLabel',{'-50','-25','0','25','50','75'})

xlabel('\Delta GPP (PgC)')
ylabel('Latitude')

set(gca,'FontSize',24)

l = legend([pH2,pH3]);
set(l,'box','off','FontSize',20)

set(gcf, 'PaperPositionMode', 'auto');
fname=strcat('./figures/diffGPP_byLatitude');

if saveFigures==1
    print(fig1,'-dpdf','-bestfit', fname)
end

%%
filename=strcat('./data/MODISlandCoverTypes.mat');
tmp=load(filename);
PFTs=tmp.dataOut; 
PFTs=PFTs(:,3:end);

% delete DNF
PFTs(PFTs==3)=NaN;
% Merge DBF and MF
PFTs(PFTs==5)=4;
% merge shrubs
PFTs(PFTs==7)=6;
% merge savana
PFTs(PFTs==9)=8;
% delete urban
PFTs(PFTs==13)=NaN;
% merge crops
PFTs(PFTs==14)=12;
% delete snow and barren
[PFTs_mapGrid, lon_b, lat_b] = regrid(lonlat,PFTs');
indxBarren_gridRS = PFTs_mapGrid>=15;



%% plot maps
%% RS
ensembleDiffAbsRS = nanmean(diffAbsRS,2);
ensembleDiffAbsRS(ensembleDiffAbsRS <-10)=-10;

Yupper = prctile(ensembleDiffAbsRS(:),90);
ensembleDiffAbsRS(ensembleDiffAbsRS>Yupper)=Yupper;

%% map the RS absolute difference
dataToMap=ensembleDiffAbsRS';
[ensembleDiffAbsRS_mapData, lon_b, lat_b] = regrid(lonlat,dataToMap);

% set a max 145 gC
ensembleDiffAbsRS_mapData(1,1)=mapMax;

mapData = ensembleDiffAbsRS_mapData;
mapData(indxBarren_gridRS)=NaN;
mapData_RS=mapData;
[lonx,latx]=meshgrid(-179.5:0.5:180,89.5:-0.5:-90);
fig1=figure;
m_proj('Equidistant cylindrical','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
m_coast('linewidth',1,'color',[0.2 0.2 0.2]); % coast line settings
hold on;

m_pcolor(lonx,latx,mapData); % draw your map here
shading INTERP;   % can be flat or INTERP
% shading flat;
colorscheme=parula(16);
colorscheme(1,:)=[0.8,0.8,0.8];

set(gca,'XTickLabel','')
colormap(colorscheme);
h=colorbar; cbfreeze(h)

set(h,'YTick',[-10,0:20:140],'ticklabels',{'<0','0','20','40','60','80','100','120','140'})
ylabel(h,'\Delta GPP (gC m^{-2} yr^{-1})','FontSize',14)
plot([-4,4],[0,0],'k--')

m_grid('box','on','tickdir','in','xticklabels',['';'']);

set(gcf, 'PaperPositionMode', 'auto');
fname=strcat('./figures/fig3c_diffAbs_rs_map_v6');


if saveFigures==1
    eval(printStatement)
    exportgraphics(fig1,strcat(fname,'.pdf'),'ContentType','vector')
    exportgraphics(fig1,strcat(fname,'.eps'),'ContentType','vector')

end








