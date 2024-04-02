% this script will...
% load and plot Sland ∆GPP relationship
% derive and plot inferred constraint

clearvars
close all

saveFigures =0; % set to 1 to save figures

addpath('./functions')
printStatement='print(fig1,''-dpdf'', fname)';

falpha=0.05;  % CI

% load GCP data
GCP_Sland = readtable('./data/GCPTerrestrialSinkSheet.xls');
startYear = 1982;

y = GCP_Sland.GCPresidualSink';
indX1 = find(GCP_Sland.Year==startYear);
indX2 = find(GCP_Sland.Year==2016);
tmp = y(indX1:indX2);
GCPbudgetBaseline(1,1:2) = [nansum(tmp), nansum(tmp)*0.185];

% 0.185 is derived from Friedlingstein et al. 2022 Table 8
% land sink estimates of 135+-25 (mean +- 1std. dev)
% 25/135 = 0.185

% load the beta values for TRENDY models, MODIS, MPI
load('./dataIntermediates/RValues_v6.mat')

%% 1. plot original R-GPP vs the Sink
modelz = [];
xdata=[];
ydata=[];
for ii=1:size(RValues.TRENDY,1)
    cModel = RValues.TRENDYnames{ii};
    modelz{ii}=cModel;

    if strcmp(cModel,'LPX')
        cModel2='LPX_Bern';
        tmp = GCP_Sland.(cModel2)(indX1:indX2);
        xdata(ii) = nansum(tmp);
    elseif strcmp(cModel,'VEGAS')
        xdata(ii) = NaN;
    elseif strcmp(cModel,'CLASS-CTEM')
        cModel2='CLASS_CTEM';
        tmp = GCP_Sland.(cModel2)(indX1:indX2);
        xdata(ii) = nansum(tmp);
    elseif strcmp(cModel,'LPJ-GUESS')
        cModel2='LPJ_GUESS';
        tmp = GCP_Sland.(cModel2)(indX1:indX2);
        xdata(ii) = nansum(tmp);
    elseif strcmp(cModel,'ORCHIDEE-MICT')
        cModel2='ORCHIDEE_MICT';
        tmp = GCP_Sland.(cModel2)(indX1:indX2);
        xdata(ii) = nansum(tmp);
     else
        cModel2 = cModel;
        tmp = GCP_Sland.(cModel2)(indX1:indX2);
        xdata(ii) = nansum(tmp);
    end

    % get Beta GPP
    ydata(ii) = RValues.TRENDY(ii,1);
end
xdata(14)=51.51; % from processing of v6 data files 

xdata2=xdata;
ydata2=ydata;

fig1=figure;

hold on
plot(xdata,ydata,'.','MarkerSize',10)
[r, prob] =corrcoef(xdata,ydata,'rows','complete');
xlabel(strcat('Cumulative S_{LAND} (',num2str(startYear),':2016)'))
ylabel('\beta_R^{GPP}')
set(gca,'FontSize',18)

% get the linear regression
x = xdata;
y = ydata;
indX=isnan(x);
x(indX)=[];y(indX)=[];
[poly, S] = polyfit(x,y,1);
slightExt = 0.05;
x2 = (min(x)-slightExt):0.01:(max(x)+slightExt);
xfit = x2;
alpha = 0.05;	% 95% CI Significance level
[Y,DELTA] = polyconf(poly,xfit,S,'alpha',alpha,'simopt','off','predopt','curve');
    
% convert CI to SE
DELTA = DELTA*2/3.92;
% plot the regression to the TRENDY models
hconf = plot(xfit,Y+DELTA,'r--');
plot(xfit,Y-DELTA,'r--')
shadedplot(xfit, Y+DELTA, Y-DELTA, [1 0.7 0.7],'none');

plot(xfit,Y+DELTA,'r--');
plot(xfit,Y-DELTA,'r--');

yfit = polyval(poly,xfit);
[Y2,DELTA2] = polyval(poly,xfit,S);
plot(xfit,Y2+DELTA2,'r--', 'color',[0.2 0.2 0.2]);
plot(xfit,Y2-DELTA2,'r--', 'color',[0.2 0.2 0.2])

yfit = polyval(poly,xfit);
p(2)=plot(x2,yfit,'r-','LineWidth',1);

% then add the actual data
plot(x,y,'k.','MarkerSize',20)

z_names = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O'};
nudgeX = 1*ones(size(xdata2));
nudgeY = 0.01*ones(size(ydata2));
nudgeY([2,3,14])=0.035;
nudgeY([1,4,5,7,8,10])=-0.02;
nudgeY(11)=0.035;
xdata2=xdata2+nudgeX;
ydata2=ydata2+nudgeY;
text(xdata2,ydata2,z_names,'FontSize',14)

% add the vertical constraint
p(1)=plot([GCPbudgetBaseline(1,1) GCPbudgetBaseline(1,1)],[0 1],'k--');
upper = GCPbudgetBaseline(1,1)+GCPbudgetBaseline(1,2);
lower = GCPbudgetBaseline(1,1)-GCPbudgetBaseline(1,2);
p(1)=plot([upper upper],[0 1],'--','color',[0.8 0.8 0.8]);
p(1)=plot([lower lower],[0 1],'--','color',[0.8 0.8 0.8]);

r_text1 = strcat('r = ',{' '},num2str(round(r(1,2),2)));
p_text2 = strcat('p = ',{' '},num2str(round(prob(1,2),2)));
r_text2 = strcat('p < ',{' '},num2str(0.01));
r_text = {r_text1{:};r_text2{:}};
text(120.25,0.32,r_text,'FontSize',18)

disp(strcat('r(Beta_{GPP}, Sland): ', num2str(r(1,2))))
disp(strcat('p(Beta_{GPP}, Sland): ', num2str(prob(1,2))))

xlim([min(xdata)-10 max(xdata)+10])
ylim([0.2 1.0])
if saveFigures==1
    fname='./figures/betaGPPVsCumSland';
    eval(printStatement)
end


%% 2. Generate the emergent constraint from the betaGPP ~ Sland relationship
nboot = 1000; %100000 used in final
nmodels = size(xdata,2);

Sland_distribution= normrnd(GCPbudgetBaseline(1,1),GCPbudgetBaseline(1,2),[1,nboot]);
nboot_obs   = size(Sland_distribution, 2);    % size of bootstrap sample of observations

xdata2=xdata;
ydata2=ydata;
nmodels2 = sum(~isnan(xdata2(1,:)));

% regression estimate of ECS PDF by bootstrap
ECSboot = zeros(nboot, 1);
for j=1:nboot
    ib      = randsample(nmodels2, nmodels2, 1);
    % UNC: uncertainty in relationship between models and Sland
    [B, BINT, R, RINT, STATS] = regress(ydata2(ib)', [ones(nmodels2, 1) xdata2(ib)']);
    sigma      = sqrt(STATS(4));

    % predicted ECS (with normal errors)
    % this bootstraps the observed Sland
    ECSboot(j) = B(1)+B(2)*Sland_distribution(1, randsample(nboot_obs, 1))+ sigma*randn(1);
end
ECSboot    = sort(ECSboot);


%% 3. Plot the inferred constraint implied by the updated RS estimates and the EC
fig1=figure;
hold on
bandwidth=0.01;
x2= 0:bandwidth:1; % the range of the x-axis
y1=[ydata, RValues.Jena(1,1), RValues.MODIS(1,1)];
y1(y1==0)=NaN;

y2 = normpdf(x2,nanmean(y1),nanstd(y1));
h1=histogram(y1,10,'Normalization','pdf','Visible','off');
b1 = bar(h1.BinEdges(1:end-1)+h1.BinWidth/2, max(y2*bandwidth)*h1.Values/max(h1.Values));
set(b1,'FaceColor',[0.8 0.8 0.8])
% plot the unconstrained PDF
p1=plot(x2,y2*bandwidth,'k.','LineWidth',2);
h = area(x2,y2*bandwidth);
h(1).FaceColor = [0.4 0.4 0.4];
h(1).FaceAlpha = 0.4;
h(1).EdgeAlpha= 0;

% now add the distribution without the RS data
y1b=ydata;
y1b(y1b==0)=NaN;

ydata_no_rs = y1b;

y2 = normpdf(x2,nanmean(ydata_no_rs),nanstd(ydata_no_rs));

% plot the unconstrained PDF
p1=plot(x2,y2*bandwidth,'k--','LineWidth',2);
h = area(x2,y2*bandwidth);
h(1).FaceColor = [0.4 0.4 0.4];
h(1).FaceAlpha = 0.4;

% now add the distribution of constrained GPP % change
jointMean = (nanmean(ECSboot)); 

varS1 = nanstd(ECSboot)^2;
jointStd2 = sqrt(varS1);

y3 = normpdf(x2,jointMean,jointStd2);
p1=plot(x2,y3*bandwidth,'-','LineWidth',2);

colorConstrained = [0.9290 0.5940 0.1250];
set(p1,'color',colorConstrained)
h = area(x2,y3*bandwidth);
h(1).FaceColor = colorConstrained;
h(1).FaceAlpha = 0.4;

text(0.3, 0.03,{'Constrained';' distribution'},'color',colorConstrained,'FontSize',18)

xlabel('\beta_R^{GPP} (%_{GPP} %_{CO2}^{-1})')
ylabel('Probability density')

set(gca,'FontSize',18)
xlim([0 1])
ylim([0 0.035])

if saveFigures==1
    fname=strcat('./figures/constrainedBeta');
    eval(printStatement)    
end

