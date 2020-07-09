function [DatesUnique,im_monthly,SheetNames,ColumnFields,jm,DateVector] = BBergDataPrecondition(filename,datesfilename,startdate,enddate)
%% Precondition Bloomberg Data
load(filename)
DatesTable = readtable(datesfilename);

Dates = table2cell(DatesTable);
FullColumnPos = find(m == max(m));

DatesDT = datetime(Dates,'InputFormat','MM/dd/yyyy','format','dd/MM/yyyy');

DatesUnique = flipud(unique(DatesDT));
DatesUnique(isnat(DatesUnique)) = [];

% Vector with first day of each month from 29/02/1996 to 11/11/2018
datedata = datenum('29/02/1996','dd/mm/yyyy'):datenum('11/11/2018','dd/mm/yyyy');
dvec   = datevec(datedata);
dvec2 = dvec(find(dvec(:,3) == 1),:);
dvecDT = datetime(dvec2,'format','dd/MM/yyyy');

[~,~,b] = intersect(DatesUnique,dvecDT);
dvecDT2 = dvecDT;
dvecDT2(b) = [];
DatesUnique = sort([DatesUnique;dvecDT2]);
clear b


im = nan(length(DatesUnique),size(MatData,2),size(MatData,3));
for i = 1:size(MatData,3)
%     disp(i)
    [s,a,b] = intersect(DatesDT(:,i),flipud(DatesUnique));
    im(b,:,i) = MatData(a,:,i);
    clear s a b
end

im = flip(im,1); % Time from going downward 1996 ---> 2018

for i = 1:size(MatData,3)
%     disp(i)
    if size(find(~isnan(im(:,1,i)),1,'first'),1) > 0
        PX_FirstNaNPos(i) = find(~isnan(im(:,1,i)),1,'first');
        PX_LastNaNPos(i) = find(~isnan(im(:,1,i)),1,'last');
    else
        PX_FirstNaNPos(i) = NaN;
        PX_LastNaNPos(i) = NaN;
    end
end
PX_FirstNaNPos = PX_FirstNaNPos';
PX_LastNaNPos = PX_LastNaNPos';

for i = 1:size(MatData,3)
%     disp(i)
    if ~isnan(PX_FirstNaNPos(i)) == 1
        im(PX_FirstNaNPos(i):PX_LastNaNPos(i),:,i) = fill_nans(im(PX_FirstNaNPos(i):PX_LastNaNPos(i),:,i));
    end
end

%% Choose first of each month from im Data
[DateVector, a, b] = intersect(DatesUnique,dvecDT);
im_monthly = im(a,:,:);  %Array contains data at first of every month.

%% Remove stocks with no PX values
NoPXStocksToRemove = find(isnan(PX_FirstNaNPos))';
SheetNames(NoPXStocksToRemove) = [];
im_monthly(:,:,NoPXStocksToRemove) = [];
im(:,:,NoPXStocksToRemove) = [];


%% Load Sector Information
jm = table(SectorNames);

    % Remove Utilities and Financials from im and sector data
    jm(NoPXStocksToRemove,:) = [];
    SectorNames(NoPXStocksToRemove) = [];
    SectorsToRemove = find(ismember(SectorNames,{'Utilities','Financials'}));
    jm(SectorsToRemove,:) = [];
    SectorNames(SectorsToRemove) = [];
    im(:,:,SectorsToRemove) = [];
    im_monthly(:,:,SectorsToRemove) = [];
    SheetNames(SectorsToRemove) = [];
    
    B2MPos = find(strcmp(ColumnFields,'PX_TO_BOOK_RATIO')); 
    %Convert M2B to B2M
    im(:,B2MPos,:) = 1./im(:,B2MPos,:);
    im_monthly(:,B2MPos,:) = 1./im_monthly(:,B2MPos,:);
    
%% Delete clone stocks
B = cellfun(@(x) x(1:end-10), SheetNames, 'un', 0);
C = unique(B);
for i = fliplr(1:length(C))
    D = find(strcmp(B,C(i)));
    if size(D,1) > 1
        D(1) = [];
        SheetNames(D) = [];
        im(:,:,D) = [];
        im_monthly(:,:,D) = [];
        jm(D,:) = [];
    end
end


%% Re-Flip im and im_monthly
im = flip(im,1);
im_monthly = flip(im_monthly,1);

    clearvars -except DatesUnique im_monthly SheetNames ColumnFields jm DateVector Dates DatesDT
    
end



