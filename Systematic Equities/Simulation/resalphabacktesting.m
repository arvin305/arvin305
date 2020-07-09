%% Monthly Residual Alpha Calculating Function
function [TSorted,ResAlpha,Alpha_zscored,Beta,bint,r,rint,stats] = resalphabacktesting(TableIn)
[m,~] = size(TableIn);

%% Beta Factors
BPPos = ReturnColumnNumber(TableIn,'B2M_SC');
MCPos = ReturnColumnNumber(TableIn,'MV_SC');
BetaPos = ReturnColumnNumber(TableIn,'Beta_SC');
DYPos = ReturnColumnNumber(TableIn,'DY_SC');

BP = TableIn.B2M_SC;
MC = TableIn.MV_SC;
Beta = TableIn.Beta_SC;
DY = TableIn.DY_SC;

%% Alpha Factors
ROAPos = ReturnColumnNumber(TableIn,'RoA_SC');
EYPos = ReturnColumnNumber(TableIn,'EY_SC');
PI6MPos = ReturnColumnNumber(TableIn,'PI6M_SC');

ROA = TableIn.RoA_SC;
EY = TableIn.EY_SC;
PI6M = TableIn.PI6M_SC;

%% Sector Risk Factors Dummy Variables
SectorPos = ReturnColumnNumber(TableIn,'Sector');
SectorVars = unique(TableIn.Sector);
for i = 1:length(SectorVars)
    str = SectorVars{i} ;
    idx = isstrprop(str,'upper') ;
    if ~exist(str(idx),'var')
        SecVarNames{i} = str(idx);
    else
        SecVarNames{i} = [str(idx),num2str(2)];
    end
end
for i = 1:length(SectorVars)
    IdenticalSVN = strcmp(SecVarNames,SecVarNames{i});
    nSVN = sum(IdenticalSVN);
    if nSVN > 1
        for j = 1:nSVN
            IdentSVNind = find(IdenticalSVN == 1);
            feval(@() evalin('caller',['SecVarNames{IdentSVNind(j)} = [SecVarNames{IdentSVNind(j)},num2str(j)];']));
        end
    end
end

for i = 1:length(SecVarNames)
    SecColumn = zeros(m,1);
    for j = 1:m
        if strcmp(SectorVars{i},TableIn.Sector{j}) == 1
            SecColumn(j) = 1;
        end
    end
    feval(@() evalin('caller',[SecVarNames{i},' = SecColumn;']));
end
SecVarStr = [SecVarNames{1}];
for i = 1:length(SectorVars)
    if i>1
        SecVarStr = [SecVarStr,',',SecVarNames{i}];
    end
end

%% Z-Scoring of factors by Sector 
INDs = {'BP','MC','ROA','EY','DY','PI6M'};
for j = 1:length(INDs)
    feval(@() evalin('caller',['',INDs{j},'bySector = zeros(m,length(SecVarNames)) ;']));
    for i = 1:length(SecVarNames)
        feval(@() evalin('caller',['',INDs{j},'bySector(:,i) = ',INDs{j},'.*',SecVarNames{i},' ; ']));
        if feval(@() evalin('caller',[' sum(',INDs{j},'bySector(:,i)) > 0 ;'])) 
            feval(@() evalin('caller',['AvgBySector = mean(',INDs{j},'bySector(',INDs{j},'bySector(:,i) > 0,i)) ; ']));
            feval(@() evalin('caller',['StdDevBySector =  std(',INDs{j},'bySector(',INDs{j},'bySector(:,i) > 0,i)) ; ']));
            if StdDevBySector ~= 0
                feval(@() evalin('caller',['',INDs{j},'bySector(',INDs{j},'bySector(:,i) > 0,i) = (',INDs{j},'bySector(',INDs{j},'bySector(:,i) > 0,i) - AvgBySector) / StdDevBySector ; ']));
            end
        end
    end
    feval(@() evalin('caller',['',INDs{j},'_zscored = sum(',INDs{j},'bySector,2); ']));
end

%% Clean up Outliers
for j = 1:length(INDs) + 1
    if j  < length(INDs) + 1
        feval(@() evalin('caller',['y = sort(',INDs{j},'_zscored);']));
    else
        feval(@() evalin('caller','y = sort(Beta);'));
    end    
    Q(1) = median(y(find(y<median(y)))); Q(2) = median(y);  Q(3) = median(y(find(y>median(y)))); IQR = Q(3)-Q(1);
    LowerOutlier = Q(1) - 1.5*IQR;
    UpperOutlier = Q(3) + 1.5*IQR;
    if j  < length(INDs) + 1
        feval(@() evalin('caller',['',INDs{j},'_zscored_cleaned = ',INDs{j},'_zscored ;']));
        feval(@() evalin('caller',['',INDs{j},'_zscored_cleaned(',INDs{j},'_zscored_cleaned < LowerOutlier) = LowerOutlier ;']));
        feval(@() evalin('caller',['',INDs{j},'_zscored_cleaned(',INDs{j},'_zscored_cleaned > UpperOutlier) = UpperOutlier ;']));
    else
        feval(@() evalin('caller','Beta_cleaned = Beta;'));
        feval(@() evalin('caller','Beta_cleaned(Beta_cleaned  < LowerOutlier) = LowerOutlier;'));
        feval(@() evalin('caller','Beta_cleaned(Beta_cleaned  > UpperOutlier) = UpperOutlier;'));
    end
end

%% Alpha Calculation
% Alpha = RankColumn(ROA_zscored) + RankColumn(EY_zscored) + RankColumn(PI6M_zscored);
Alpha = 0.4*(ROA_zscored) + 0.6*(EY_zscored + PI6M_zscored);
Alpha_zscored = (Alpha - mean(Alpha)) / std(Alpha);

%% Compute Coefficients through multiple linear regression
feval(@() evalin('caller',['BetaIndVars = [BP_zscored,MC_zscored,Beta,DY_zscored,',SecVarStr,'];'])); %ones(size(Alpha)),
[Beta,bint,r,rint,stats] = regress(Alpha_zscored,BetaIndVars);

%% Residual Alpha computation and Table Ranking
ResAlpha = Alpha_zscored - BetaIndVars*Beta;
ResAlphaRanked = RankColumn(ResAlpha);
ResAlphaIndex = [];
for i = 1:m
    ResAlphaIndex = [ResAlphaIndex;find(ResAlphaRanked == i)];
end
ResAlphaIndex = ResAlphaIndex';
TSorted = TableIn(ResAlphaIndex,:);

end

% % Plot Beta Regression
% scatter3(BP,MC,Alpha,'filled')
% hold on
% BPfit = linspace(min(BP),max(BP),20);
% MCfit = linspace(min(MC),max(MC),20);
% [BPFIT,MCFIT] = meshgrid(BPfit,MCfit);
% ALPHAFIT = Beta(1) + Beta(2)*BPFIT + Beta(3)*MCFIT;
% mesh(BPFIT,MCFIT,ALPHAFIT)
% xlabel('BP')
% ylabel('MC')
% zlabel('ALPHA')
% hold off