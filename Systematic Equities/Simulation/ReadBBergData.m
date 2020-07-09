%% Use Code to Convert Bloomberg Data from Excel Sheet to 3D Array
clc; clear;
exl = actxserver('excel.application');
exlWkbk = exl.Workbooks;
exlFile = exlWkbk.Open([pwd,'\RIY_2018_11_11_001-800.xlsx']);
exlSheets = exlFile.Sheets;
% Get number of open Workbooks
numSheets = exlSheets.Count;
% Query name of each open Workbook
SheetNames = cell(1,numSheets);
for i = 1:numSheets
    SheetNames{i} = exlSheets.Item(i).Name;
end

disp('file open!')

% Data = zeros(NumberofRows,20,numSheets);
Data = struct();
for i = 1:numSheets
    disp(i);
    exlSheet1 = exlSheets.Item(SheetNames{i});
    m(i) = exlSheet1.Columns.End(4).row;       % Find the end of the column
    n(i) = exlSheet1.Columns.End(2).column;    % Find the end of the column
    rngObjFull = exlSheet1.Range(['A1:',xlscol(n(i)),num2str(m(i))]);
    evalin('base',['Data.Stock',num2str(i),' = rngObjFull.Value;']);  %Save Raw data as a Cell
    evalin('base',['Data.Stock',num2str(i),'{2,1} = ''11/11/2018'' ;']);  %Save Raw data as a Cell
end

MatData = nan(max(m)-1,20,numSheets);
DateData = cell(max(m)-1,numSheets);
for i = 1:numSheets
    disp(i)
    % Extract Numeric Data
    evalin('base',['Out = Data.Stock',num2str(i),'(2:end,2:n);']);  %Save Raw data as a Cell
    idn = cellfun(@isnumeric,Out); % identify numeric values.
    out = nan(size(Out,1),20);
    out(idn) = [Out{idn}];         % allocate numeric values.
    MatData(1:size(out,1),:,i) = out;
    
    %Extract Date Strings  
    DateOut = {};
    evalin('base',['DateOut = Data.Stock',num2str(i),'(2:m(i),1);']);  %Save Raw data as a Cell
    DateData(1:m(i)-1,i) = DateOut;        
end
DateDataTable = cell2table(DateData);
writetable(DateDataTable,'DateData1.xlsx')

% Data = permute(Data,[1 3 2]);
% Data(:,:,9)

% Use code below to convert / into space if there are any in Sector names
% A = SheetNames  
% B = RIYSectorNames(:,1)
% B = strrep(B,'/',' ')

Quit(exl);
% delete(exl);