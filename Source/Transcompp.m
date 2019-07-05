%{ 
This file is part of TRANSCOMPP
Copyright (C) 2019  N. Suhas Jagannathan

TRANSCOMPP is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

TRANSCOMPP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with TRANSCOMPP.  If not, see <https://www.gnu.org/licenses/>.
%}

function varargout = Transcompp(varargin)

%Initialization code
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @transCompp_OpeningFcn, ...
    'gui_OutputFcn',  @transCompp_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before transCompp is made visible.
function transCompp_OpeningFcn(hObject, eventdata, handles, varargin)
setAllPanelsVisibility(hObject, handles,'off');
set(handles.loadDataPanel,'visible','on');
handles = resetAllData(hObject, eventdata, handles);

% --- Outputs from this function are returned to the command line.
function varargout = transCompp_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --------------------------------------------------------------------
function fileMenu_Callback(hObject, eventdata, handles)


function handles = resetAllData(hObject, eventdata, handles)
handles.output = hObject;
handles.replicateFiles = {};
handles.replicateList = {};
handles.replicateData = {};
handles.scData = {};
handles.popFracDataFromSC = [];
handles.scatter = [];
handles.popFracFile = {};
handles.popFracData = [];

handles.simParams = {};
handles.simParams.errorStringVal = 1;
handles.simParams.numIterOpt = 50;
handles.simParams.numSamples = 500;
handles.simParams.resample = 0;
handles.simParams.bdConstr = [];
handles.simParams.transConstr = [];
handles.simParams.solveForBD = 0;
handles.simParams.lts_frac = 0.8;
handles.simParams.scResampleSize = 100;
handles.simParams.scStateDef = {};

handles.inputData.dataToUse = [];
handles.inputData.timepoints = [];
handles.inputData.numEnt = [];

handles.results.transMat = [];
handles.results.bd = [];

inequalityString = {'<';'<=';'=';'>';'>='};
set(handles.inequalityDropdown,'String',inequalityString);
set(handles.LoadFracFileButton,'Enable','off');
set(handles.customErrorFileChooseButton,'enable','inactive');
set(handles.replicateFileInfoText,'String',' ');
set(handles.bdStateVal,'String',' ');
set(handles.bdLBVal,'String',' ');
set(handles.bdUBval,'String',' ');
set(handles.BDConstraintsText,'String',' ');
set(handles.transConstrState1Val,'String',' ');
set(handles.transConstrState2Val,'String',' ');
set(handles.transConstrLBval,'String',' ');
set(handles.transConstrUBval,'String',' ');
set(handles.transConstraintsText,'String',' ');
set(handles.currStateText,'String',' ');
set(handles.transConstraintsText,'String',' ');
set(handles.defineStateThresholdVal,'String',' ');
set(handles.definedStatesText,'String',' ');
set(handles.TrimFractionVal,'String',' ');
set(handles.numIterBestFitVal,'String',' ');
set(handles.numResamplesVal,'String',' ');
set(handles.transMatRowForHistVal,'String',' ');
set(handles.transMatColForHistVal,'String',' ');
set(handles.resampleCheckBox,'Value',0);

set(handles.selectReplicateDropdown,'String',' ');
set(handles.xAxisDropdown,'String',' ');
set(handles.yAxisDropdown,'String',' ');
set(handles.addStatesTable,'Data',[]);

set(handles.bestFitTransMatTable,'data',cell(2,2));
cla(handles.transMatHistAxes,'reset');

handles.currentState = '';

set(handles.summaryDataTypeVal,'string',' ');
set(handles.summaryNumStatesVal,'string',' ');
set(handles.summaryNumTimepointsVal,'string',' ');
set(handles.summaryDataTable,'Data',[]);
set(handles.summarySolveForBDVal,'String','No');
set(handles.summaryNumBDConstrVal,'String','0');
set(handles.summaryNumTransProbConstrVal,'String','0');
set(handles.summaryResamplingVal,'String','Off');

set(handles.popFracDataFilenameText, 'string',' ');
cla(handles.axesSCData,'reset');
set(handles.xAxisMinVal,'string',' ');
set(handles.xAxisMaxVal,'string',' ');
set(handles.yAxisMinVal,'string',' ');
set(handles.yAxisMaxVal,'string',' ');
set(handles.xAxisLogScaleCheckbox,'value',0);
set(handles.yAxisLogScaleCheckbox,'value',0);
set(handles.solveBDCheckbox,'value',0);
guidata(hObject, handles);

% --------------------------------------------------------------------
function newAnalysisMenuItem_Callback(hObject, eventdata, handles)
opts.Interpreter = 'tex';
opts.Default = 'Cancel';
answer = questdlg('\fontsize{9} This will erase all loaded data. Are you sure you want to continue?',...
    'New Analysis',...
    'Continue','Cancel',opts);
switch answer
    case 'Continue'
        handles = resetAllData(hObject, eventdata, handles);
        setAllPanelsVisibility(hObject, handles,'off');
        set(handles.loadDataPanel,'visible','on');
    case 'Cancel'
        return;
end
guidata(hObject,handles);

% --------------------------------------------------------------------
function loadAnalysisMenuItem_Callback(hObject, eventdata, handles)
global h;
handles_orig = handles;
try
[filename, foldername] = uigetfile({'*.mat'}, 'Select an analysis file');
eval(['load ',foldername,filename,'']);
handles.replicateFiles = analysis.replicateFiles;
handles.replicateList = analysis.replicateList;
handles.replicateData = analysis.replicateData;
handles.scData = analysis.scData;
handles.popFracDataFromSC = analysis.popFracDataFromSC;
handles.scatter = analysis.scatter;
handles.popFracFile = analysis.popFracFile;
handles.popFracData = analysis.popFracData;
handles.simParams = analysis.simParams;
handles.inputData = analysis.inputData;
handles.results = analysis.results;

if ~isempty(handles.replicateData)
    set(handles.selectReplicateDropdown,'string',handles.replicateList);
    set(handles.selectReplicateDropdown,'Value',1);
    selectReplicateDropdown_Callback(hObject, eventdata, handles); 
end

if handles.simParams.solveForBD
    set(handles.solveBDCheckbox,'Enable','on');
end

if handles.simParams.resample    
    set(handles.resampleCheckBox,'Enable','on');
end

switch handles.simParams.errorStringVal
    case 1
        set(handles.SSERadioButton,'Enable','on');
    case 2
        set(handles.L1radioButton,'Enable','on');
    case 3
        set(handles.TLSradioButton,'Enable','on');
        set(handles.TrimFractionVal,'string',num2str(handles.simParams.lts_frac));
    case 4
        set(handles.customErrorRadioButton,'Enable','on');
end

set(handles.numIterBestFitVal,'string',num2str(handles.simParams.numIterOpt));
set(handles.numResamplesVal,'string',num2str(handles.simParams.numSamples));

if ~isempty(handles.simParams.bdConstr)
    b = handles.simParams.bdConstr;
    bd_str = {};
    for i = 1:size(b,1)
        state = b(i,1);
        bd_lb = b(i,2);
        bd_ub = b(i,3);
        strToAdd = ['Added BD constraint for state ',num2str(state),' : ','LB = ',num2str(bd_lb),' UB = ',num2str(bd_ub)];
        bd_str = [bd_str; {strToAdd}];
    end   
    set(handles.BDConstraintsText,'String',bd_str);
end

if ~isempty(handles.simParams.transConstr)
    t = handles.simParams.transConstr;
    trans_str = {};
    for i = 1:size(t,1)
        state1 = t(i,1);
        state2 = t(i,2);
        trans_lb = t(i,3);
        trans_ub = t(i,4);
        strToAdd = ['TransProb b/w states ',num2str(state1),',',num2str(state2),' : ','LB = ',num2str(trans_lb),' UB = ',num2str(trans_ub)];
        trans_str = [trans_str; {strToAdd}];
    end
    set(handles.transConstraintsText,'String',trans_str);
end

if ~isempty(handles.simParams.scStateDef)
    stateText = {};
    numStates = size(handles.simParams.scStateDef,1);
    for i = 1:numStates
        currState = handles.simParams.scStateDef{i,1};
        state_str = stateString(currState);
        stateText = [stateText; {state_str}];
    end
    set(handles.definedStatesText,'string',stateText);
end
catch
    errordlg('Please check the input file','Load failed','modal');
    handles = handles_orig;
    return;
end
msgbox('Successfully loaded analysis file','Load Success','modal');
guidata(hObject,handles);
h = handles;


% --------------------------------------------------------------------
function saveAnalysisMenuItem_Callback(hObject, eventdata, handles)
analysis.replicateFiles = handles.replicateFiles;
analysis.replicateList = handles.replicateList;
analysis.replicateData = handles.replicateData;
analysis.scData = handles.scData;
analysis.popFracDataFromSC = handles.popFracDataFromSC;
analysis.scatter = handles.scatter;
analysis.popFracFile = handles.popFracFile;
analysis.popFracData = handles.popFracData;

analysis.simParams = handles.simParams;
analysis.inputData = handles.inputData;
analysis.results = handles.results;
uisave('analysis');

% --------------------------------------------------------------------
function exitMenuItem_Callback(hObject, eventdata, handles)
close(handles.figure1);


% --- Executes during object creation, after setting all properties.
function replicateFileInfoText_CreateFcn(hObject, eventdata, handles)
set(hObject,'HorizontalAlignment','left');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addReplicateButton.
function addReplicateButton_Callback(hObject, eventdata, handles)
global h;
[filename, foldername] = uigetfile({'*.xlsx','*.xls'}, 'Select a replicate data file');
replicateNum = length(handles.replicateFiles)+1;
try
    [~,sheets] = xlsfinfo([foldername,filename]);
    if ~isempty(sheets)
        timepoints = [];
        for i = 1:length(sheets)
            t = str2double(sheets{i});
            if isnan(t)
                errordlg('Please ensure that the sheet names are all numeric values corresponding to the timepoints','Sheet Name Error','modal');
                handles = clearAllSCData(hObject, eventdata, handles);
                h = handles;
                return;
            end    
            
            if mod(t,1) > 0
                errordlg('Sheet names can only be integers','Non-integer Timepoint Error','modal');
                set(handles.popFracDataFilenameText,'String',' ');
                h = handles;
                return;
            end
            
            [~,~,tmp] = xlsread([foldername,filename],sheets{i});
            if size(tmp,1) < 101
                errordlg('Too few cells. At least 100 cells needed for single cell data','SC data error','modal');
                handles = clearAllSCData(hObject, eventdata, handles);
                h = handles;
                return;
            end
            
            
            cols = tmp(1,:);
            if ~iscellstr(cols)
                errordlg('Please ensure that the first row in each worksheet has the column name','Column Name Error','modal');
                handles = clearAllSCData(hObject, eventdata, handles);
                h = handles;
                return;
            end
            
            if i == 1
                base_cols = cols;
            end
            
            if ~isequal(base_cols,cols)
                errordlg('Please ensure that the column names are in the same order in each sheet (timepoint)','Column Name Mismatch Error','modal');
                handles = clearAllSCData(hObject, eventdata, handles);
                h = handles;
                return;
            end
                
            handles.replicateData{replicateNum,i} = tmp;
            timepoints = [timepoints t];
        end
    end
   
    handles.inputData.timepoints = timepoints;
    handles.replicateList = [handles.replicateList;['Replicate ',num2str(replicateNum)]];

    set(handles.selectReplicateDropdown,'string',handles.replicateList);
    selectReplicateDropdown_Callback(hObject, eventdata, handles);
    handles.replicateFiles = [handles.replicateFiles;filename];
    replicateListString = get(handles.replicateFileInfoText,'string');
    strToAdd = ['Replicate ',num2str(replicateNum),' : ',filename];
    if isequal(replicateListString,' ')
        replicateListString = strToAdd;
    else
        replicateListString = [replicateListString;{strToAdd}];
    end
    set(handles.replicateFileInfoText,'string',replicateListString);
    handles.popFracData = [];
    handles.popFracFile = {};
    set(handles.popFracDataFilenameText,'string',' ');  
    set(handles.summaryDataTypeVal,'string','Single-Cell Data');
    set(handles.summaryNumTimepointsVal,'string',num2str(length(timepoints)));
    set(handles.summaryDataTable,'data',[]);
catch
end
h = handles;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function addReplicateButton_CreateFcn(hObject, eventdata, handles)

% --- Executes on selection change in popupmenuReplicateList.
function popupmenuReplicateList_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenuReplicateList_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TrimFractionVal_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function TrimFractionVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function numIterBestFitVal_Callback(hObject, eventdata, handles)
try
    n = str2double(get(hObject,'String'));
    if ~isnan(n)
        handles.simParams.numIterOpt = n;
    end
catch
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function numIterBestFitVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in resampleCheckBox.
function resampleCheckBox_Callback(hObject, eventdata, handles)
try
    n = get(hObject,'Value');
    if ~isnan(n)
        handles.simParams.resample = n;
        if n == 1
            set(handles.summaryResamplingVal,'String','On');
        elseif n == 0
            set(handles.summaryResamplingVal,'String','Off');
        end
    end
catch
end
guidata(hObject,handles);

function numResamplesVal_Callback(hObject, eventdata, handles)
try
    n = str2double(get(hObject,'String'));
    if ~isnan(n)
        handles.simParams.numSamples = n;
    end
catch
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function numResamplesVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in defaultSimParamsButton.
function defaultSimParamsButton_Callback(hObject, eventdata, handles)
set(handles.numIterBestFitVal,'String','50');
set(handles.numResamplesVal,'String','500');
set(handles.resampleCheckBox,'Value',0);
guidata(hObject,handles);

function transMatRowForHistVal_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function transMatRowForHistVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function transMatColForHistVal_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function transMatColForHistVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in errorTermButtonGroup.
function errorTermButtonGroup_SelectionChangeFcn(hObject, eventdata, handles)
errorString = get(eventdata.NewValue,'String');
trimFrac = 0.8;
errorStringVal = 1;
if strcmpi(errorString,'Least Squares')
    errorStringVal = 1;
    set(handles.customErrorFileChooseButton,'enable','inactive');
elseif strcmpi(errorString,'L1 Norm')
    errorStringVal = 2;
    set(handles.customErrorFileChooseButton,'enable','inactive');
elseif strcmpi(errorString,'Trimmed least squares')
    errorStringVal = 3;
    try
        trimFrac = str2double(get(handles.TrimFractionVal,'String'));
        set(handles.customErrorFileChooseButton,'enable','inactive');
    catch
    end
elseif strcmpi(errorString,'Custom error file')
    errorStringVal = 4;
    set(handles.customErrorFileChooseButton,'enable','on');
    disp('Please enter yor own error files');
end
handles.simParams.errorStringVal = errorStringVal;
handles.simParams.lts_frac = trimFrac;
guidata(hObject,handles);

% --- Executes on button press in LoadFracFileButton.
function LoadFracFileButton_Callback(hObject, eventdata, handles)
global h;
[filename, foldername] = uigetfile({'*.xlsx','*.xls'}, 'Select a replicate data file');
set(handles.popFracDataFilenameText,'String',[foldername,filename]);
[~,sheets] = xlsfinfo([foldername,filename]);
if ~isempty(sheets)
    popfrac = []; timepoints = [];
    for i = 1:length(sheets)
        t = str2double(sheets{i});
        if isnan(t)
            errordlg('Please ensure that the sheet names are all numeric values corresponding to the timepoints','Sheet Name Error','modal');
            set(handles.popFracDataFilenameText,'String',' ');
            h = handles;
            return;
        end
        if mod(t,1) > 0
            errordlg('Sheet names can only be integers','Non-integer Timepoint Error','modal');
            set(handles.popFracDataFilenameText,'String',' ');
            h = handles;
            return;
        end
             
        [~,~,tmp] = xlsread([foldername,filename],sheets{i});

        try
            popfrac = [popfrac,tmp];
        catch
            errordlg('Every worksheet should have the same number of elements','Inconsistent size error','modal')
        end
        numEnt = size(tmp,2);
        timepoints = [timepoints t];       
    end
end

tmp = []; row = 1; sort = 1;
merged_popfrac = [];
for i = 1:size(popfrac,1)
    if i > 1  && mean(isnan(cell2mat(popfrac(i,:)))) == 1
        row = 1;
        sort = sort + 1;
        continue;
    end
    tmp(row,:,sort) = cell2mat(popfrac(i,:));
    merged_popfrac = [merged_popfrac; cell2mat(popfrac(i,:))];
    row = row  + 1;
end

handles.popFracData = tmp;
handles.inputData.dataToUse = tmp;
handles.inputData.timepoints = timepoints;
handles.inputData.numEnt = numEnt;
handles = clearAllSCData(hObject, eventdata, handles);
set(handles.summaryDataTypeVal,'string','Population fraction data');
set(handles.summaryNumStatesVal,'string',num2str(numEnt));
set(handles.summaryNumTimepointsVal,'string',num2str(length(timepoints)));
set(handles.summaryDataTable,'data',merged_popfrac);
guidata(hObject,handles);

function handles = clearAllSCData(hObject, eventdata, handles)
handles.scData = {};
set(handles.selectReplicateDropdown,'string','Select Replicate');
set(handles.xAxisDropdown,'string',' ');
set(handles.yAxisDropdown,'string',' ');
set(handles.xAxisDropdown,'string','X axis');
set(handles.yAxisDropdown,'string','Y axis');
set(handles.xAxisMinVal,'string',' ');
set(handles.xAxisMaxVal,'string',' ');
set(handles.yAxisMinVal,'string',' ');
set(handles.yAxisMaxVal,'string',' ');
set(handles.defineStateThresholdVal,'string',' ');
set(handles.replicateFileInfoText,'string',' ');
set(handles.definedStatesText,'string',' ');
set(handles.addStatesTable,'Data',cell(2,2));
cla(handles.axesSCData,'reset');
set(handles.selectReplicateDropdown,'string',' ');
handles.replicateFiles = {};
handles.replicateList = {}; 
handles.replicateData = {};
handles.popFracDataFromSC = [];
guidata(hObject,handles);

function popFracDataFilenameText_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popFracDataFilenameText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in runSimButton.
function runSimButton_Callback(hObject, eventdata, handles)
global h;
h = handles;
numEntities = handles.inputData.numEnt;
timepoints = handles.inputData.timepoints;
numTriesOpt = handles.simParams.numIterOpt;
trans_constr = handles.simParams.transConstr;
solveForBD = handles.simParams.solveForBD;
bd_constr =  handles.simParams.bdConstr;
errorTerm = handles.simParams.errorStringVal;
lts_frac = str2double(get(handles.TrimFractionVal,'String'));
numSim = handles.simParams.numSamples;
bootstrapSize = handles.simParams.scResampleSize;
to_resample = handles.simParams.resample
scStateDef = [];

if handles.simParams.solveForBD
    bd_vector = 1;
else
    if isempty(bd_constr)
        bd_vector = 1;
    else
        bd_vector = ones(numEntities,1);
        for i = 1:size(bd_constr,1)
            state = bd_constr(i,1);
            lb = bd_constr(i,2);
            ub = bd_constr(i,3);
            bd_vector(state) = mean(lb,ub);
        end
    end
end

if ~to_resample
    solveType = 1;
    if ~isempty(handles.popFracDataFromSC)
      TimeSeriesData = handles.popFracDataFromSC; 
    elseif ~isempty(handles.popFracData)
      TimeSeriesData = handles.popFracData; 
    else
        errordlg ({'Sorry there is no population fraction data available to carry out this operation.';
                'If you have loaded single cell data, please make sure you define states and generate population fractions'},'Input data error','modal');
        return;
    end
else
    if ~isempty(handles.scData)
        solveType = 2;
        TimeSeriesData = handles.scData;   
        scStateDef = handles.simParams.scStateDef;
    elseif ~isempty(handles.popFracData)
        solveType = 3;
        TimeSeriesData = handles.popFracData;
    elseif ~isempty(handles.popFracDataFromSC)
        solveType = 3;
        TimeSeriesData = handles.popFracDataFromSC;
    else
        errordlg('Please check your input data type and define states if necessary','Input state error','modal');
        return;
    end
end
solveType;
TimeSeriesData;
bd_vector;
[transMat,bd_solved] = runSimulation(numEntities,TimeSeriesData,timepoints,numTriesOpt,trans_constr,solveForBD,bd_vector,bd_constr,errorTerm,lts_frac,solveType,numSim,bootstrapSize,scStateDef);
transMat;
bd_solved;

if to_resample
    mean_TM = mean(transMat);
    std_TM = std(transMat);
    str_tmp = cell(size(mean_TM));
    for i = 1:length(mean_TM)
        str_tmp{i} = [num2str(mean_TM(i),3),' +/- ',num2str(std_TM(i),3)];
    end
    str_tmp = reshape(str_tmp,handles.inputData.numEnt,handles.inputData.numEnt)';
else
    str_tmp = transMat;
end
set(handles.bestFitTransMatTable,'Data',str_tmp);
handles.results.transMat = transMat;
handles.results.bd = bd_solved;
msgbox('Simulation complete!','Success','modal');
guidata(hObject,handles);


% --- Executes on button press in plotTransMatHistButton.
function plotTransMatHistButton_Callback(hObject, eventdata, handles)
try
    r = str2double(get(handles.transMatRowForHistVal,'String'));
    c = str2double(get(handles.transMatColForHistVal,'String'));
    i = (r-1)*handles.inputData.numEnt+c;
    axes(handles.transMatHistAxes);
    hist(handles.results.transMat(:,i),25);
catch
end

function numStatesDefineVal_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function numStatesDefineVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu9.
function popupmenu9_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenu9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function transConstraintsText_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function transConstraintsText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function transConstrState1Val_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function transConstrState1Val_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function transConstrState2Val_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function transConstrState2Val_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function transConstrLBval_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function transConstrLBval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function transConstrUBval_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function transConstrUBval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in transConstraintAddButton.
function transConstraintAddButton_Callback(hObject, eventdata, handles)
try
    if isempty(handles.inputData.numEnt) || handles.inputData.numEnt == 0
        errordlg('If using Single cell data, Please define states first','Trans Constr error','modal');
        return;
    end
    trans_constr = handles.simParams.transConstr;
    state1 = str2double(get(handles.transConstrState1Val,'String'));
    state2 = str2double(get(handles.transConstrState2Val,'String'));
    trans_lb = str2double(get(handles.transConstrLBval,'String'));
    trans_ub = str2double(get(handles.transConstrUBval,'String'));
    trans_constr = [trans_constr; state1 state2 trans_lb trans_ub];
    
    if mod(state1,1) > 0 || mod(state2,1) > 0 || state1 > handles.inputData.numEnt || state2 > handles.inputData.numEnt
        errordlg('Incorrect State numbers specified','State error','modal');
        return;
    end
    
    if trans_lb < 0 || trans_lb > 1 || trans_ub < 0 || trans_ub > 1
        errordlg('Transition rates can only be between 0 and 1','Transition rate error','modal');
        return;
    end
    
    handles.simParams.transConstr = trans_constr;
    trans_str = get(handles.transConstraintsText,'String');
    strToAdd = ['TransProb b/w states ',num2str(state1),',',num2str(state2),' : ','LB = ',num2str(trans_lb),' UB = ',num2str(trans_ub)];
    trans_str = [trans_str; {strToAdd}];
    set(handles.transConstraintsText,'String',trans_str); 
    set(handles.summaryNumTransProbConstrVal,'string',num2str(size(trans_constr,1)));
catch
    disp('Error in entered Transition rate constraints');
end
guidata(hObject,handles)

function BDConstraintsText_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function BDConstraintsText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bdStateVal_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function bdStateVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bdLBVal_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function bdLBVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function dataMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function loadDataMenu_Callback(hObject, eventdata, handles)
setAllPanelsVisibility(hObject,handles,'off');
set(handles.loadDataPanel,'visible','on');
guidata(hObject,handles);


% --------------------------------------------------------------------
function simulationMenu_Callback(hObject, eventdata, handles)

% --- Executes on button press in scDataRadioButton.
function scDataRadioButton_Callback(hObject, eventdata, handles)
scData = get(hObject,'Value');
if scData
    set(handles.addReplicateButton,'Enable','on');
    set(handles.LoadFracFileButton,'Enable','off');
    set(handles.replicateFileInfoText,'Enable','on');
end
guidata(hObject,handles);

% --- Executes on button press in popFracDataRadioButton.
function popFracDataRadioButton_Callback(hObject, eventdata, handles)
fracData = get(hObject,'Value');
if fracData
    set(handles.addReplicateButton,'Enable','off');
    set(handles.LoadFracFileButton,'Enable','on');
    set(handles.replicateFileInfoText,'Enable','off');
end
guidata(hObject,handles);

% --- Executes on button press in clearDataLoadPageButton.
function clearDataLoadPageButton_Callback(hObject, eventdata, handles)
handles.popFracFile = {};
handles.scatter = [];
handles.popFracData = [];
handles.inputData.dataToUse = [];
handles.inputData.timepoints = [];
handles.inputData.numEnt = [];
set(handles.popFracDataFilenameText,'String',' ');
handles = clearAllSCData(hObject, eventdata, handles);
set(handles.bestFitTransMatTable,'data',cell(2,2));
cla(handles.transMatHistAxes,'reset');
handles.replicateFiles = {};
set(handles.selectReplicateDropdown,'string',' ');
set(handles.summaryDataTypeVal,'string',' ');
set(handles.summaryNumStatesVal,'string',' ');
set(handles.summaryNumTimepointsVal,'string',' ');
set(handles.summaryDataTable,'data',[]);
guidata(hObject,handles);

function bdUBval_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function bdUBval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function scData = getSCDataFromReplicateData(handles)
rData = handles.replicateData;
for i = 1:size(rData,1)
    for j = 1:size(rData,2)
        tmp = rData{i,j};
        if ischar(tmp{1,1})
            tmp(1,:) = [];
        end
        m =  cell2mat(tmp);
        [~,p] = pca(m);
        m = [m p];
        rData{i,j} = m;
    end
end
scData = rData;

% --- Executes on button press in genPopFracButton.
function genPopFracButton_Callback(hObject, eventdata, handles)
global h;
rData = getSCDataFromReplicateData(handles);
handles.scData = rData;
popFracData = [];
for i = 1:size(rData,1)
    pop_frac_row = [];
    for j = 1:size(rData,2)
        tmp = rData{i,j}';
        a = handles.simParams.scStateDef;
        if isempty(handles.simParams.scStateDef)
            errordlg('Please enter state definitions before generating population fractions','State Def Error','modal');
            return;
        end

        pop_frac = getPopFrac(handles.simParams.scStateDef,tmp);
        pop_frac_row = [pop_frac_row pop_frac];
    end
    popFracData(i,:) = pop_frac_row;
end
handles.popFracDataFromSC = popFracData;
handles.inputData.numEnt = size(handles.simParams.scStateDef,1);
currText = get(handles.definedStatesText,'string');
if isequal(currText, ' ')
    currText = {};
end
strToAdd = 'Population Fractions Generated!';
currText = [currText; {strToAdd}];
set(handles.definedStatesText,'string',currText);
set(handles.summaryNumStatesVal,'string',num2str(handles.inputData.numEnt));
set(handles.summaryDataTable,'data',popFracData);
h = handles;
guidata(hObject,handles);

function definedStatesText_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function definedStatesText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in addStateDefButton.
function addStateDefButton_Callback(hObject, eventdata, handles)
try
    global h;
    currState = handles.currentState;
    tmp = readStateDef(hObject, eventdata, handles);
    if sum(tmp{1,1} ~= 0)
        currState = [currState; tmp];
    end
    state_str = stateString(currState);
    numStates = size(handles.simParams.scStateDef,1);
    handles.simParams.scStateDef{numStates+1,1} = currState;
    currText = get(handles.definedStatesText,'string');
    if isequal(currText, ' ')
        currText = {};
    end
    currText = [currText; {state_str}];
    set(handles.definedStatesText,'string',currText);
    handles.currentState = '';
    h = handles;
catch
end
set(handles.currStateText,'String',' ');
guidata(hObject,handles);

function tmp = readStateDef(hObject, eventdata, handles)
try
    tmp = {};
    coeff = get(handles.addStatesTable,'Data'); coeff = cell2mat(coeff(:,2))';
    ineq = get(handles.inequalityDropdown,'Value');
    thres = str2double(get(handles.defineStateThresholdVal,'string'));
    tmp{1,1} = coeff; tmp{1,2} = ineq; tmp{1,3} = thres; tmp{1,4} = '';
catch
end

% --- Executes on button press in clearBDButton.
function clearBDButton_Callback(hObject, eventdata, handles)
handles.simParams.bdConstr = [];
set(handles.BDConstraintsText,'String',' ');
set(handles.bdStateVal,'String',' ');
set(handles.bdLBVal,'String',' ');
set(handles.bdUBval,'String',' ');
set(handles.solveBDCheckbox,'value',0)
set(handles.summarySolveForBDVal,'string','No');
set(handles.summaryNumBDConstrVal,'string','0');
handles.simParams.solveForBD = 0;
guidata(hObject,handles);

% --- Executes on button press in clearTransConstraintButton.
function clearTransConstraintButton_Callback(hObject, eventdata, handles)
handles.simParams.transConstr = [];
set(handles.transConstrState1Val,'String',' ');
set(handles.transConstrState2Val,'String',' ');
set(handles.transConstrLBval,'String',' ');
set(handles.transConstrUBval,'String',' ');
set(handles.transConstraintsText,'String',' ');
set(handles.summaryNumTransProbConstrVal,'string','0');
guidata(hObject,handles);

% --- Executes on selection change in selectReplicateDropdown.
function selectReplicateDropdown_Callback(hObject, eventdata, handles)
val = get(handles.selectReplicateDropdown,'Value');
data = handles.replicateData{val};
cols = data(1,:)';
l = min(length(cols),3);
pc_legends = {'PC1';'PC2';'PC3'};
cols = [cols;pc_legends(1:l)];
set(handles.xAxisDropdown,'string',cols);
set(handles.yAxisDropdown,'string',cols);
stateDef = cell(length(cols),2);
stateDef(:,1) = cols;
stateDef(:,2) = num2cell(zeros(length(cols),1));
set(handles.addStatesTable,'Data',stateDef);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function selectReplicateDropdown_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in xAxisDropdown.
function xAxisDropdown_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function xAxisDropdown_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in yAxisDropdown.
function yAxisDropdown_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function yAxisDropdown_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotSCDataButton.
function plotSCDataButton_Callback(hObject, eventdata, handles)
pc_legends = {'PC1','PC2','PC3'};
xcol = get(handles.xAxisDropdown,'value');
ycol = get(handles.yAxisDropdown,'value');
rNum = get(handles.selectReplicateDropdown,'value');
data = handles.replicateData{rNum};
c = size(data,2);
m =  cell2mat(data(2:end,:));
l = 1;
labels = data(1,:);

if xcol > c || ycol > c
    [~,p] = pca(m);
    l = min(size(p,2),3);
    m = [m p];
    labels = [labels,pc_legends(1:l)];
end

axes(handles.axesSCData);
scatter(m(:,xcol),m(:,ycol),'.');


xlabel(labels(xcol));
ylabel(labels(ycol));
handles.scatter = [m(:,xcol) m(:,ycol)];
guidata(hObject,handles);

function xAxisMinVal_Callback(hObject, eventdata, handles)
x_lim = xlim;
try
    val = str2double(get(hObject,'String'));
    xlim([val x_lim(2)]);
catch
    xlim([min(handles.scatter(:,1)) max(handles.scatter(:,1))]);
end

% --- Executes during object creation, after setting all properties.
function xAxisMinVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xAxisMaxVal_Callback(hObject, eventdata, handles)
x_lim = xlim;
try
    val = str2double(get(hObject,'String'));
    xlim([x_lim(1) val]);
catch
    xlim([min(handles.scatter(:,1)) max(handles.scatter(:,1))]);
end

% --- Executes during object creation, after setting all properties.
function xAxisMaxVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function yAxisMinVal_Callback(hObject, eventdata, handles)
y_lim = ylim;
try
    val = str2double(get(hObject,'String'));
    ylim([val y_lim(2)]);
catch
    ylim([min(handles.scatter(:,2)) max(handles.scatter(:,2))]);
end

% --- Executes during object creation, after setting all properties.
function yAxisMinVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function yAxisMaxVal_Callback(hObject, eventdata, handles)
y_lim = ylim;
try
    val = str2double(get(hObject,'String'));
    ylim([y_lim(1) val]);
catch
    ylim([min(handles.scatter(:,2)) max(handles.scatter(:,2))]);
end

% --- Executes during object creation, after setting all properties.
function yAxisMaxVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in xAxisLogScaleCheckbox.
function xAxisLogScaleCheckbox_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
if val == 1
    set(gca,'Xscale','log');
else
    set(gca,'Xscale','linear');
end

% --- Executes on button press in yAxisLogScaleCheckbox.
function yAxisLogScaleCheckbox_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
if val == 1
    set(gca,'Yscale','log');
else
    set(gca,'Yscale','linear');
end

% --- Executes on button press in clearStateDefButton.
function clearStateDefButton_Callback(hObject, eventdata, handles)
handles.simParams.scStateDef = {};
set(handles.definedStatesText,'String',' ');
set(handles.summaryNumStatesVal,'String',' ');
set(handles.summaryDataTable,'data',[]);
guidata(hObject,handles);

% --- Executes on button press in solveBDCheckbox.
function solveBDCheckbox_Callback(hObject, eventdata, handles)
n = get(hObject,'Value');
handles.simParams.solveForBD = n;
if n == 0
    set(handles.summarySolveForBDVal,'string','No');
elseif n == 1
    set(handles.summarySolveForBDVal,'string','Yes');
end
guidata(hObject,handles);

function setAllPanelsVisibility(hObject, handles,visibility)
set(handles.loadDataPanel,'visible',visibility);
set(handles.modelConstraintsPanel,'visible',visibility);
set(handles.runSimulationPanel,'visible',visibility);
set(handles.visualizeSCPanel,'visible',visibility);
set(handles.viewResultsPanel,'visible',visibility);
set(handles.defineStatesPanel,'visible',visibility);
guidata(hObject,handles);

% --------------------------------------------------------------------
function panelsMenuItem_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function showPanelsMenuItem_Callback(hObject, eventdata, handles)
setAllPanelsVisibility(hObject, handles,'on');
guidata(hObject,handles);

% --------------------------------------------------------------------
function hidePanelsMenuItem_Callback(hObject, eventdata, handles)
setAllPanelsVisibility(hObject, handles,'off');
guidata(hObject,handles);

% --------------------------------------------------------------------
function modelMenuItem_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function modelConstraintsMenuItem_Callback(hObject, eventdata, handles)
setAllPanelsVisibility(hObject, handles,'off');
set(handles.modelConstraintsPanel,'visible','on');
guidata(hObject,handles);

% --------------------------------------------------------------------
function defineStatesMenuItem_Callback(hObject, eventdata, handles)
setAllPanelsVisibility(hObject, handles,'off');
set(handles.defineStatesPanel,'visible','on');
guidata(hObject,handles);

% --------------------------------------------------------------------
function visualizeSCDataMenuItem_Callback(hObject, eventdata, handles)
setAllPanelsVisibility(hObject, handles,'off');
set(handles.visualizeSCPanel,'visible','on');
guidata(hObject,handles);

% --------------------------------------------------------------------
function runSimulationMenuItem_Callback(hObject, eventdata, handles)
setAllPanelsVisibility(hObject, handles,'off');
set(handles.runSimulationPanel,'visible','on');
guidata(hObject,handles)

% --------------------------------------------------------------------
function viewResultsMenuItem_Callback(hObject, eventdata, handles)
setAllPanelsVisibility(hObject, handles,'off');
set(handles.viewResultsPanel,'visible','on');
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function axesSCData_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in bdConstraintAddButton.
function bdConstraintAddButton_Callback(hObject, eventdata, handles)
try
    if isempty(handles.inputData.numEnt) || handles.inputData.numEnt == 0
        errordlg('If using Single cell data, Please define states first','BD error','modal');
        return;
    end
    bd_constr = handles.simParams.bdConstr;
    state = str2double(get(handles.bdStateVal,'String'));
    bd_lb = str2double(get(handles.bdLBVal,'String'));
    bd_ub = str2double(get(handles.bdUBval,'String'));
    bd_constr = [bd_constr; state bd_lb bd_ub];
    if mod(state,1) > 0 || state > handles.inputData.numEnt
        errordlg('Invalid state parameter entered','State number error','modal');
        return;
    end 
    handles.simParams.bdConstr = bd_constr;
    bd_str = get(handles.BDConstraintsText,'String');
    strToAdd = ['BD constraint for state ',num2str(state),' : ','LB = ',num2str(bd_lb),' UB = ',num2str(bd_ub)];
    bd_str = [bd_str; {strToAdd}];
    set(handles.BDConstraintsText,'String',bd_str);
    set(handles.summaryNumBDConstrVal,'string',num2str(size(bd_constr,1)));
catch
    disp('Error in entered BD constraints');
end
guidata(hObject,handles)


% --------------------------------------------------------------------
function errorTermButtonGroup_ButtonDownFcn(hObject, eventdata, handles)

function runSimButton_ButtonDownFcn(hObject, eventdata, handles)

function defineStateThresholdVal_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function defineStateThresholdVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in inequalityDropdown.
function inequalityDropdown_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function inequalityDropdown_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
     
% --- Executes on button press in andButton.
function andButton_Callback(hObject, eventdata, handles)
global h;
currState = handles.currentState;
tmp = readStateDef(hObject, eventdata, handles);
tmp{1,4} = [tmp{1,4},'and'];
handles.currentState = [currState; tmp];
currState_str = getCurrentStateString(handles.currentState);
set(handles.currStateText,'String',currState_str);
coeff = get(handles.addStatesTable,'Data'); coeff(:,2) = num2cell(0*cell2mat(coeff(:,2)));
set(handles.addStatesTable,'Data',coeff);
h = handles;
guidata(hObject,handles);

% --- Executes on button press in orButton.
function orButton_Callback(hObject, eventdata, handles)
currState = handles.currentState;
tmp = readStateDef(hObject, eventdata, handles);
tmp{1,4} = [tmp{1,4},'or'];
handles.currentState = [currState; tmp];
currState_str = getCurrentStateString(handles.currentState);
set(handles.currStateText,'String',currState_str);
coeff = get(handles.addStatesTable,'Data'); coeff(:,2) = num2cell(0*cell2mat(coeff(:,2)));
set(handles.addStatesTable,'Data',coeff);
guidata(hObject,handles);

% --- Executes on button press in openParButton.
function openParButton_Callback(hObject, eventdata, handles)
global h;
currState = handles.currentState;
if isequal(currState,'')
    warndlg('Adding State: ');
end
l = size(currState,1);
currState{l,4} = [currState{l,4},'('];
handles.currentState = currState;
currState_str = getCurrentStateString(handles.currentState);
set(handles.currStateText,'String',currState_str);
h = handles;
guidata(hObject,handles);

% --- Executes on button press in closeParButton.
function closeParButton_Callback(hObject, eventdata, handles)
currState = handles.currentState;
tmp = readStateDef(hObject, eventdata, handles);
tmp{1,4} = [tmp{1,4},')'];
handles.currentState = [currState; tmp];
currState_str = getCurrentStateString(handles.currentState);
set(handles.currStateText,'String',currState_str);
coeff = get(handles.addStatesTable,'Data'); coeff(:,2) = num2cell(0*cell2mat(coeff(:,2)));
set(handles.addStatesTable,'Data',coeff);
guidata(hObject,handles);

function state_str = stateString(state)
l = size(state,1);
all_ineq = {'<','<=','=','>','>='};
state_str = 'Added state : ';
for i = 1:l
    tmp = state(i,:);
    strToAdd = ['[',num2str(tmp{1,1}),']*scData ',all_ineq{tmp{1,2}},' ',num2str(tmp{1,3}),' ',tmp{1,4},' '];
    state_str = [state_str strToAdd];    
end

function currStateText_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function currStateText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in cancelCurrentStateButton.
function cancelCurrentStateButton_Callback(hObject, eventdata, handles)
set(handles.currStateText,'String',' ');
handles.currentState = '';
guidata(hObject,handles);

function currState_str = getCurrentStateString(currentState)
ineq_vals = {'<','<=','==','>','>='};
    l = size(currentState,1);
    currState_str = '';
    for i = 1:l
        coeff = currentState{i,1};
        ineq = currentState{i,2};
        thres = currentState{i,3};
        chain = currentState{i,4};
        to_add = ['[',num2str(coeff),']*scData ',' ',ineq_vals{ineq},' ',num2str(thres),' ',chain];
        currState_str = [currState_str,to_add];
    end

% --- Executes on button press in customErrorFileChooseButton.
function customErrorFileChooseButton_Callback(hObject, eventdata, handles)
[errorFileName, errorFileDir] = uigetfile('*.m');
errorFile = [errorFileDir,errorFileName];
try
copyfile(errorFile,'customError.m');
catch
    errordlg('Please choose a valid matlab function file','CustomFileError','modal');
    waitforbuttonpress
    customErrorFileChooseButton_Callback(hObject, eventdata, handles)
end
msgbox('Custom Error file added successfully','Success','modal');
guidata(hObject,handles);


% --------------------------------------------------------------------
function exportResultsMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to exportResultsMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global h;
h = handles;
results = handles.results;
uisave('results');
