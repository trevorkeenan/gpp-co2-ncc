function [A, lon,lat]= regrid(lonlatIn,valueIn)

% converrt the p model style lat lon list to a grid

lonlat=lonlatIn;

% set the dimensions / lonlat of the entire grid
lat=-89.75:0.5:89.75;
lon=-179.75:0.5:179.75;
m=length(lat);
n=length(lon);

lonlat2=[];
for ii = 1:m
    latz=lat(ii)*ones(n,1);
    lonz=lon';
    tmp=[lonz,latz];
    lonlat2=vertcat(lonlat2,tmp);
end
clear lonz latz ii tmp 

% find location of land in whole grid
indXClimate=ismember(lonlat2,lonlat,'rows');
[lon,lat]=meshgrid(-179.75:0.5:179.75,89.75:-0.5:-89.75);

%% Regrid the data
tmp=nan(length(indXClimate),1);
indX=valueIn==0;
neighbour=valueIn';
neighbour=vertcat(neighbour(2:end),0);
valueIn(indX)=neighbour(indX);
tmp(indXClimate)=valueIn;
A=reshape(tmp,n,m);

 mapdata=A;
        mapdata(mapdata==0)=NaN;
        mapdata=flipud(mapdata'); % a 2 D` matrix
A=mapdata;
end