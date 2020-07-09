%% Backtester
clear; clc; close all;

%% Load Bloomberg Data Stocks Data
filename = 'HISTCOMP_RIY_2018_11_11.mat';
datesfilename = 'Dates.xlsx';
startdate = '29/02/1996';
enddate = '11/11/2018';

% [DatesUnique,im,SheetNames,ColumnFields,jm,DateVector] = BBergDataPrecondition(filename,datesfilename,startdate,enddate);
% If preconditioned data exists load mat file instead of above
load('PreconditionedData.mat')

%% Find variable Positions in Column Fields
RoAPos = find(strcmp(ColumnFields,'RETURN_ON_ASSET')); 
B2MPos = find(strcmp(ColumnFields,'PX_TO_BOOK_RATIO')); 
MVPos = find(strcmp(ColumnFields,'HISTORICAL_MARKET_CAP')); 
BetaPos = find(strcmp(ColumnFields,'BETA_RAW_OVERRIDABLE')); 
EYPos = find(strcmp(ColumnFields,'EBIT_EV_YIELD')); 
DYPos = find(strcmp(ColumnFields,'AVERAGE_DIVIDEND_YIELD')); 
ACPos = find(strcmp(ColumnFields,'PX_LAST'));

%% Parameters
Period = 12;

%% Variables
StocksZ = im;
AC = permute(StocksZ(:,ACPos,:),[1 3 2]); AC_flipped = flipud(AC);
DY = permute(StocksZ(:,DYPos,:),[1 3 2]); %DY_flipped = flipud(DY);
B2M = permute(StocksZ(:,B2MPos,:),[1 3 2]); %B2M_flipped = flipud(B2M);
MV = permute(StocksZ(:,MVPos,:),[1 3 2]); %MV_flipped = flipud(MV);
RoA = permute(StocksZ(:,RoAPos,:),[1 3 2]); %RoA_flipped = flipud(RoA);
EY = permute(StocksZ(:,EYPos,:),[1 3 2]); %EY_flipped = flipud(EY);
Beta = permute(StocksZ(:,BetaPos,:),[1 3 2]); %Beta_flipped = flipud(Beta);

[imx, imy, imz] = size(im);

%% Load S&P 500 Data
filename = 'S&P500.xlsx';
SP500 = readtable(filename);
SP500.Properties.VariableNames(1) = {'Date'};

formatIn = 'mmm dd yyyy';
% SP500Dates = string(SP500.Date);
% SP500Dates = datetime(datevec(SP500Dates,formatIn),'format','dd/MM/yyyy');
SP500Dates = datetime(datevec(SP500.Date,formatIn),'format','dd/MM/yyyy');

[~, ~, b] = intersect(SP500Dates,DateVector);
DateIndexVector = (1:size(DateVector));
DateIndexVector(b) = [];
clear s a b;

SP500Dates2 = sort([SP500Dates;DateVector(DateIndexVector)]);
[~,loc] = ismember(SP500Dates,SP500Dates2);

SP500Data = nan(size(SP500Dates2));
SP500Data(loc) = SP500.Price;
SP500Data = fill_nans(SP500Data);

[s, a, b] = intersect(SP500Dates2,DateVector);
SP500Dates = DateVector;
SP500Data = SP500Data(a);


%% Calculate Forward Returns
ACRestriction = isnan(AC); % if 1 then NaN for Adjusted Close Value (Do not trade at NaN value)
ACRestriction_flipped = flipud(ACRestriction);
n = [1,3,6,12];

for k = 1:length(n)
    FR_Flipped = zeros(imx-n(k),imz);
    for i = 1:imz
        for j = 1:imx
            if j > imx-n(k)
                FR_Flipped(j,i) = NaN;                
            else
                if isnan(AC_flipped(j,i)) == 1 || isnan(AC_flipped(j + n(k), i )) == 1
                    FR_Flipped(j,i) = NaN;
                else if isnan(AC_flipped(j + n(k), i )) == 1
                        FR_Flipped(j,i) = 0;
                    else
                        FR_Flipped(j,i) = (AC_flipped(j + n(k) , i) - AC_flipped(j, i)) ./ AC_flipped(j,i);
                    end
                end
            end
        end
    end
    evalin('base',['FR',num2str(n(k)),' = flipud(FR_Flipped);']);
end

for i = imy+1:imy+4
    evalin('base',['im(:,i,:) = FR',num2str(n(i-imy)),';']);
    evalin('base',['FR',num2str(n(i-imy)),'_Flipped = flipud(FR',num2str(n(i-imy)),');']);
end

PI6M = circshift(FR6,[-1 0]);
PI6M_Flipped = flipud(PI6M);


%% Define Which Stocks Can be traded When using Boolean matrix:
TradeRestrictionBoolean = ACRestriction; % TradeBoolean calculated for tradeable stocks as per AC
% Debugging: B2M_Test = B2M(:,[36 148]');AC_Test = AC(:,[36 148]');TradeRestrictionBoolean = TradeRestrictionBoolean(:,[36 148]');[TradeRestrictionBoolean_Out,B2M_Test_Out_flipped,B2M_Test_Out] = IndicatorDataPreconditioning(TradeRestrictionBoolean,B2M_Test,AC_Test);[TradeRestrictionBoolean_Out,B2M_Test_Out,AC_Test]
IND = {'B2M','MV','RoA','Beta','EY','DY','PI6M'};

for i = 1:length(IND)
    evalin('base',['[TradeRestrictionBoolean,',IND{i},'_flipped,',IND{i},'] = IndicatorDataPreconditioning(TradeRestrictionBoolean,',IND{i},',AC);']);
end

for i = 1:length(n)
    evalin('base',['FR',num2str(n(i)),'Boolean = isnan(FR',num2str(n(i)),'); FR',num2str(n(i)),'Boolean_flipped = flipud(FR',num2str(n(i)),'Boolean);']);
end


% % Update im variable with pre-conditioned variables
% im(:,:,RoAPos) = RoA';
% im(:,:,B2MPos) = B2M';
% im(:,:,MVPos) = MV';
% im(:,:,BetaPos) = Beta';
% im(:,:,EYPos) = EY';
% im(:,:,DYPos) = DY';

%% Refresh StocksZ and MonthsZ variables with new cleaned variables and calculated FR data
StocksZ = im;
MonthsZ = permute(StocksZ,[3 2 1]);

%% Find Investable Stocks at a given month for a given investment time period
InvestmentCycle = 2:Period:imx;
SelectionCycle = 1:Period:imx;
SelectionCycle = SelectionCycle(1:length(InvestmentCycle));
TradeRestrictionBoolean_Flipped = flipud(TradeRestrictionBoolean);
evalin('base',['TradeRestrictionBoolean_Flipped_InvPeriod = TradeRestrictionBoolean_Flipped(InvestmentCycle,:) + TradeRestrictionBoolean_Flipped(SelectionCycle,:) + FR',num2str(Period),'Boolean_flipped(InvestmentCycle,:);']); % Filters the Stocks as per data availability during Selection Month, the investing month after selection month and availbility of Future Return Values for the given period
TradeRestrictionBoolean_Flipped_InvPeriod(TradeRestrictionBoolean_Flipped_InvPeriod > 1 ) = 1;
if sum(TradeRestrictionBoolean_Flipped_InvPeriod(end,:) ) == date
    InvestmentCycle(end) = [];
    SelectionCycle(end) = [];
    TradeRestrictionBoolean_Flipped_InvPeriod(end,:) = [];
end

ReturnCycle = InvestmentCycle + Period;

RawStockListBoolean_Flipped = (TradeRestrictionBoolean_Flipped_InvPeriod == 0);
ResAlpha = nan([length(SelectionCycle),max(sum(RawStockListBoolean_Flipped,2))]);
Alpha = nan([length(SelectionCycle),max(sum(RawStockListBoolean_Flipped,2))]);
Beta = nan([length(SelectionCycle),40]);


for i = 1:length(SelectionCycle)
    i
    
    RawStockList = find(RawStockListBoolean_Flipped(i,:));
    B2M_SC_flipped = B2M_flipped(SelectionCycle(i),RawStockList)'; B2M_SC = flipud(B2M_SC_flipped);
    MV_SC_flipped = MV_flipped(SelectionCycle(i),RawStockList)'; MV_SC = flipud(MV_SC_flipped);
    RoA_SC_flipped = RoA_flipped(SelectionCycle(i),RawStockList)'; RoA_SC = flipud(RoA_SC_flipped);
    Beta_SC_flipped = Beta_flipped(SelectionCycle(i),RawStockList)'; Beta_SC = flipud(Beta_SC_flipped);
    EY_SC_flipped = EY_flipped(SelectionCycle(i),RawStockList)'; EY_SC = flipud(EY_SC_flipped);
    DY_SC_flipped = DY_flipped(SelectionCycle(i),RawStockList)'; DY_SC = flipud(DY_SC_flipped);
    PI6M_SC_flipped = PI6M_flipped(SelectionCycle(i),RawStockList)'; PI6M_SC = flipud(PI6M_SC_flipped);
    
    Sector = jm.SectorNames(RawStockList);
    StockIndex = RawStockList';
    T = table(StockIndex,SheetNames(StockIndex),Sector,B2M_SC,MV_SC,RoA_SC,Beta_SC,EY_SC,DY_SC,PI6M_SC);
    if sum(RawStockListBoolean_Flipped(i,:)) > 1
        
        [TSorted,ResAlpha_i,Alpha_i,Beta_i,bint,r,rint,stats] = resalphabacktesting(T);
        Beta(i,1:size(Beta_i,1)) = Beta_i';
        Alpha(i,1:size(Alpha_i,1)) = Alpha_i';
        ResAlpha(i,1:size(ResAlpha_i,1)) = ResAlpha_i';
        
        SelectedStocks(i,:) = TSorted.StockIndex(1:10);
        SelectedStocksShorting(i,:) = TSorted.StockIndex(end-9:end);
        B2M_SelectedStocks(i,:) = TSorted.B2M_SC(1:10);
        MV_SelectedStocks(i,:) = TSorted.MV_SC(1:10);
        RoA_SelectedStocks(i,:) = TSorted.RoA_SC(1:10);
        Beta_SelectedStocks(i,:) = TSorted.Beta_SC(1:10);
        EY_SelectedStocks(i,:) = TSorted.EY_SC(1:10);
        DY_SelectedStocks(i,:) = TSorted.DY_SC(1:10);
        PI6M_SelectedStocks(i,:) = TSorted.DY_SC(1:10);

    else
        if size(T) > 0
            disp('Only one stock');
            for i=1:10
                TSorted(i,:) = T(1,:);
            end
            SelectedStocksShorting(1,:) = TSorted.StockIndex(end-10:end);
            SelectedStocks(1,:) = TSorted.StockIndex(1:10);
        else
            SelectedStocks(1,1:10) = NaN;
            SelectedStocksShorting(1,1:10) = NaN;
        end
    end
    
    size(SelectedStocks)
end

for i = 1:length(InvestmentCycle)
    if sum(isnan(SelectedStocks(i,:))) == 10
        evalin('base',['FR',num2str(Period),'_Flipped_SelectedStocks(i,1:10) = 0;']);
    else
        evalin('base',['FR',num2str(Period),'_Flipped_SelectedStocks(i,:) = FR',num2str(Period),'_Flipped(InvestmentCycle(i),SelectedStocks(i,:));']);
        evalin('base',['FR',num2str(Period),'Boolean_flipped_SelectedStocks(i,:) = FR',num2str(Period),'Boolean_flipped(InvestmentCycle(i),SelectedStocks(i,:));']);
    end
end
evalin('base',['AverageReturn = mean(FR',num2str(Period),'_Flipped_SelectedStocks,2);']);

mean(AverageReturn);

GrowthPerPeriod = AverageReturn + 1;
for i = 1:length(GrowthPerPeriod)
    TotalGrowth(i) = prod(GrowthPerPeriod(1:i));
end
% [InvestmentCycle',ReturnCycle',AverageReturn,GrowthPerPeriod,TotalGrowth']

%% Create Time Vector for Plotting
TotalGrowth = [1,TotalGrowth]'; 
GrowthPerPeriod = [1;GrowthPerPeriod]; 
TimeVector = [InvestmentCycle(1),ReturnCycle]'; 

InvestmentDateVector = DateVector(TimeVector);

%% Calculate S&P Data Returns for same Dates
SP500DataForOurTimeVector = SP500Data(TimeVector);
SP500TotalGrowth = SP500DataForOurTimeVector./SP500DataForOurTimeVector(1);
SP500GrowthPerPeriod = [1;SP500TotalGrowth(2:end)./SP500TotalGrowth(1:end-1)];

%% Sharpe Ratio
% 1-Year Return
YearSubCycles = 12/Period;
NumberOfYearsInData = floor((length(InvestmentDateVector) - 1)/YearSubCycles);
NumberOfSubCycleRows = 1 + NumberOfYearsInData*YearSubCycles;

SharpeYearVector = InvestmentDateVector((1:YearSubCycles:NumberOfSubCycleRows));
SharpeYearGrowth = TotalGrowth((1:YearSubCycles:NumberOfSubCycleRows));
SharpeYearGrowthPerPeriod = [1;SharpeYearGrowth(2:end)./SharpeYearGrowth(1:end-1)];  %Returns

Rf = (1.0234);
SharpeRatio1Year = (mean(SharpeYearGrowthPerPeriod) - Rf)/std(SharpeYearGrowthPerPeriod)

% 3-Year Return
Year3SubCycles = 36/Period;
NumberOf3YearsInData = floor((length(InvestmentDateVector) - 1)/Year3SubCycles);
NumberOf3SubCycleRows = 1 + NumberOf3YearsInData*Year3SubCycles;

Sharpe3YearVector = InvestmentDateVector((1:Year3SubCycles:NumberOf3SubCycleRows));
Sharpe3YearGrowth = TotalGrowth((1:Year3SubCycles:NumberOf3SubCycleRows));
Sharpe3YearGrowthPerPeriod = [1;Sharpe3YearGrowth(2:end)./Sharpe3YearGrowth(1:end-1)];  %Returns

Rf3 = (1.0234)^3;
SharpeRatio3Year = (mean(Sharpe3YearGrowthPerPeriod) - Rf3)/std(Sharpe3YearGrowthPerPeriod)

%% Info Ratio
SP500YearGrowth = SP500TotalGrowth((1:YearSubCycles:NumberOfSubCycleRows));
SP500YearGrowthPerPeriod = [1;SP500YearGrowth(2:end)./SP500YearGrowth(1:end-1)];  %Returns

InfoRatio = (SharpeYearGrowthPerPeriod - SP500YearGrowthPerPeriod)/std(SharpeYearGrowthPerPeriod - SP500YearGrowthPerPeriod);

disp('Annual Returns / Info Ratio')
disp([100*(SharpeYearGrowthPerPeriod-1),InfoRatio])

%% Plots
g = plot(InvestmentDateVector,TotalGrowth,InvestmentDateVector,GrowthPerPeriod,InvestmentDateVector,SP500TotalGrowth,InvestmentDateVector,SP500GrowthPerPeriod);
set(g,'LineWidth',3);
gLegend = legend('Strategy Total','Strategy Returns','S&P500 Total','S&P Returns');
set(gLegend,'location','NorthWest','FontSize',18);
set(gca, 'FontName', 'Helvetica','FontSize', 14,'Box','on', 'TickDir', 'in', 'TickLength', [.02 .02],'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'off','XColor', [ 0 0 0], 'YColor', [0 0 0],'LineWidth', 2);

StrategyAnnualReturn = (TotalGrowth(end)).^(1/NumberOfYearsInData)
SP500AnnualReturn = (SP500TotalGrowth(end)).^(1/NumberOfYearsInData)

% % Debugging
% BugStockInd = RawStockList(find(isnan(B2M_SC_flipped))); %Find Bugged Stocks from RawStockList
% InvestmentCycle(i) %Current Investment Cycle
% TradeRestrictionBoolean_Flipped(38,BugStockInd')
% B2M_flipped(:,BugStockInd')
% [(1:imy)',B2M_flipped(:,BugStockInd'),TradeRestrictionBoolean_Flipped(:,BugStockInd')]

% filename = fullfile('QI_Screener_Export_10_25_2018.csv');
% T = readtable(filename);
% T.Properties.VariableNames{1} = char('ClosingPrice');
% % Code to search for words in Table Column
% mask = ismember(YourTable{:,3}, TheWord);
% YourTable(mask,4) = YourTable(mask,5);
% CPPos = ReturnColumnNumber(T,'ClosingPrice');
% CP = T.ClosingPrice;




%% Load Data Code for Python Data
% im = readNPY('hist_factor_data 2018-11-03.pkl.npy');
% im(:,:,9:12) = [];
% im(:,:,7) = []; % X - STOCKS , Y - MONTHS, Z - VARIABLES ----> ['Book to Market', 'Market Value (USD)', 'Return on Assets', 'Beta','Earnings Yield', 'Dividend Yield', 'Adj_close', 'fret1', 'fret2', 'fret3', 'fret6'])    
% jm = readtable('SectorNames.xlsx');
% % jm([find(strcmp(jm.SectorNames,'nm'));find(strcmp(jm.SectorNames,'na'))],:) = [];
% jm = jm(1:imx,:);
