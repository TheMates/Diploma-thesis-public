function varargout = fittingGUI(varargin)
% FITTINGGUI MATLAB code for fittingGUI.fig
%      FITTINGGUI, by itself, creates a new FITTINGGUI or raises the existing
%      singleton*.
%
%      H = FITTINGGUI returns the handle to a new FITTINGGUI or the handle to
%      the existing singleton*.
%
%      FITTINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FITTINGGUI.M with the given input arguments.
%
%      FITTINGGUI('Property','Value',...) creates a new FITTINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fittingGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fittingGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fittingGUI

% Last Modified by GUIDE v2.5 08-May-2020 15:56:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fittingGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @fittingGUI_OutputFcn, ...
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


% --- Executes just before fittingGUI is made visible.
function fittingGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for fittingGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
%LOAD MY COMMON CODE
mydir  = pwd;idcs   = strfind(mydir,filesep);root_dir = mydir(1:idcs(end)-1); %one directory up
addpath([root_dir filesep 'common_code'])
handles.pathBox.String = [root_dir filesep 'Devices data'];

% UIWAIT makes fittingGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fittingGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on slider movement.
function freqStartSlider_Callback(hObject, eventdata, handles)
value = round(get(hObject,'Value'));
set(handles.freqStartBox, 'String', num2str(value));
set(handles.freqStartSlider, 'Value', value);

plotIt(handles);



% --- Executes during object creation, after setting all properties.
function freqStartSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on slider movement.
function freqStopSlider_Callback(hObject, eventdata, handles)
value = round(get(hObject,'Value'));
set(handles.freqStopBox, 'String', num2str(value));
set(handles.freqStopSlider, 'Value', value);
plotIt(handles);

% --- Executes during object creation, after setting all properties.
function freqStopSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function crossFreqSlider_Callback(hObject, eventdata, handles)
value = round(get(hObject,'Value'));
set(handles.crossFreqBox, 'String', num2str(value));
set(handles.crossFreqSlider, 'Value', value);
plotIt(handles);

% --- Executes during object creation, after setting all properties.
function crossFreqSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function crossLengthSlider_Callback(hObject, eventdata, handles)
value = round(get(hObject,'Value'));
set(handles.crossLengthBox, 'String', num2str(value));
set(handles.crossLengthSlider, 'Value', value);
plotIt(handles);

% --- Executes during object creation, after setting all properties.
function crossLengthSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function nPoles1Slider_Callback(hObject, eventdata, handles)
value = round(get(hObject,'Value'));
set(handles.nPoles1Box, 'String', num2str(value));
set(handles.nPoles1Slider, 'Value', value);
if handles.fixedNPolesCheck.Value == 1
    fixedNPoles = str2double( handles.fixedNPolesBox.String);
    set(handles.nPoles2Box, 'String', num2str(fixedNPoles-value));
    set(handles.nPoles2Slider, 'Value', fixedNPoles-value);
end
plotIt(handles);

% --- Executes during object creation, after setting all properties.
function nPoles1Slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function nPoles2Slider_Callback(hObject, eventdata, handles)
value = round(get(hObject,'Value'));
set(handles.nPoles2Box, 'String', num2str(value));
set(handles.nPoles2Slider, 'Value', value);
plotIt(handles);

% --- Executes during object creation, after setting all properties.
function nPoles2Slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function lambda1Slider_Callback(hObject, eventdata, handles)
value = get(hObject,'Value');
set(handles.lambda1Box, 'String', num2str(value));
set(handles.lambda1Slider, 'Value', value);
plotIt(handles);


% --- Executes during object creation, after setting all properties.
function lambda1Slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function lambda2Slider_Callback(hObject, eventdata, handles)
value = get(hObject,'Value');
set(handles.lambda2Box, 'String', num2str(value));
set(handles.lambda2Slider, 'Value', value);
plotIt(handles);


% --- Executes during object creation, after setting all properties.
function lambda2Slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



%%%%%%%%%%%%%%%%%%%%% BOXES %%%%%%%%%%%%%

function freqStartBox_Callback(hObject, eventdata, handles)
value = round(str2num( get(hObject,'String')));
set(handles.freqStartBox, 'String', num2str(value));
set(handles.freqStartSlider, 'Value', value);
plotIt(handles);


% --- Executes during object creation, after setting all properties.
function freqStartBox_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function freqStopBox_Callback(hObject, eventdata, handles)
value = round(str2num( get(hObject,'String')));
set(handles.freqStopBox, 'String', num2str(value));
set(handles.freqStopSlider, 'Value', value);
plotIt(handles);

% --- Executes during object creation, after setting all properties.
function freqStopBox_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function crossFreqBox_Callback(hObject, eventdata, handles)
value = round(str2num( get(hObject,'String')));
set(handles.crossFreqBox, 'String', num2str(value));
set(handles.crossFreqSlider, 'Value', value);
plotIt(handles);

% --- Executes during object creation, after setting all properties.
function crossFreqBox_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function crossLengthBox_Callback(hObject, eventdata, handles)
value = round(str2num( get(hObject,'String')));
set(handles.crossLengthBox, 'String', num2str(value));
set(handles.crossLengthSlider, 'Value', value);
plotIt(handles);

% --- Executes during object creation, after setting all properties.
function crossLengthBox_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function nPoles1Box_Callback(hObject, eventdata, handles)
value = round(str2num( get(hObject,'String')));
set(handles.nPoles1Box, 'String', num2str(value));
set(handles.nPoles1Slider, 'Value', value);
plotIt(handles);

% --- Executes during object creation, after setting all properties.
function nPoles1Box_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function nPoles2Box_Callback(hObject, eventdata, handles)
value = round(str2num( get(hObject,'String')));
set(handles.nPoles2Box, 'String', num2str(value));
set(handles.nPoles2Slider, 'Value', value);
plotIt(handles);

% --- Executes during object creation, after setting all properties.
function nPoles2Box_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function lambda1Box_Callback(hObject, eventdata, handles)
value = str2num( get(hObject,'String'));
set(handles.lambda1Box, 'String', num2str(value));
set(handles.lambda1Slider, 'Value', value);
plotIt(handles);

% --- Executes during object creation, after setting all properties.
function lambda1Box_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function lambda2Box_Callback(hObject, eventdata, handles)
value = str2num( get(hObject,'String'));
set(handles.lambda2Box, 'String', num2str(value));
set(handles.lambda2Slider, 'Value', value);
plotIt(handles);


% --- Executes during object creation, after setting all properties.
function lambda2Box_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function plotIt(handles)
list = get(handles.fileListbox,'String');
file = list{get(handles.fileListbox,'Value')};
path = get(handles.pathBox,'String');

freqStart = get(handles.freqStartSlider, 'Value');
freqStop =  get(handles.freqStopSlider, 'Value');
crossFreq =  get(handles.crossFreqSlider, 'Value');
crossLength = get(handles.crossLengthSlider, 'Value');
nPoles1 = get(handles.nPoles1Slider, 'Value');
nPoles2 = get(handles.nPoles2Slider, 'Value');
lambda1 = get(handles.lambda1Slider, 'Value');
lambda2 = get(handles.lambda2Slider, 'Value');
ITER = str2num(get(handles.iterBox, 'String'));
Fs = str2double(handles.boxFs.String);
NFIR = str2double(handles.NFIRBox.String);
nakSpline = handles.NAKSplineCheck.Value;
axes(handles.plot1);

if ~handles.viewOnlyCheck.Value
    
    load([path '\' file]);  %fr,H
    w = 2*pi*fr/Fs;
    if handles.UseCppCheck.Value
        [Bm,Am,FIR] = filterDesignCpp(w,H, freqStart, freqStop, crossFreq, crossLength, nPoles1, nPoles2,lambda1, lambda2,NFIR,ITER,Fs,nakSpline);
    else
        [Bm,Am,FIR] = filterDesignMatlab(w,H, freqStart, freqStop, crossFreq, crossLength, nPoles1, nPoles2,lambda1, lambda2,NFIR,ITER,Fs,~nakSpline);
    end
    p = drawPlot(Bm,Am,FIR,w,Fs,H);
     %vykreslení pólù
     axes(handles.plot2);
     zplane(p(abs(p)<1));
     if ~isempty(p(abs(p)>=1))
     zplane(p);
     lims = axis;
     zplane(p(abs(p)<1));
     hold on
     [hz1,hp1,ht1] = zplane(p(abs(p)>=1));
     set(findobj(hz1, 'Type', 'line'),'Color', 'r');
     set(findobj(hp1, 'Type', 'line'),'Color', 'r');
     hold off
     set(gca,'color',[1 .9 .9]);
     axis(lims);
     end
else    %view only
    load([path '\' file]);      %Bm, Am, FIR, (fr,H ?)
    if exist('fr','var') && exist('H','var')
        semilogx(fr,db(H),'g--');
        hold on;
    end
    fr = logspace(log10(freqStart), log10(freqStop),1024);fr = fr';
    w = 2*pi*fr/Fs;
    h = bankCore.parfiltfresp(Bm,Am,FIR,fr,Fs);
    semilogx(fr,db(h));
    hold off;
    
%     mydir = path;        %LOAD MY COMMON CODE
%     data_path = [mydir(1:strfind(path,'coeffs')-1) 'data\'];
%     load([data_path file]);
%     hold on
%     semilogx(fr,db(H),'g--');
%     hold off
    grid on
    axis([10 Fs/2 -50 30]);
    
end    
% --- Executes on button press in loadBtn.
function loadBtn_Callback(~, eventdata, handles)
folder_name = get(handles.pathBox,'String');
% get what is inside the folder
Infolder = dir(folder_name);
% Initialize the cell of string that will be update in the list box
MyListOfFiles = {Infolder(~[Infolder.isdir]).name};
% update the listbox with the result
set(handles.fileListbox, 'String',  MyListOfFiles);
plotIt(handles)


% --- Executes on selection change in fileListbox.
function fileListbox_Callback(hObject, eventdata, handles)
list = get(handles.fileListbox,'String');
file = list{get(handles.fileListbox,'Value')};

set(handles.exportCoeffsBox,'String',file(1:length(file)-4));

plotIt(handles)



% --- Executes during object creation, after setting all properties.
function fileListbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function iterBox_Callback(hObject, eventdata, handles)
plotIt(handles);

% --- Executes during object creation, after setting all properties.
function iterBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--- Executes on button press in exportCoeffsBox.
function exportCoeffsBox_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function exportCoeffsBox_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in exportCoeffsBtn.
function exportCoeffsBtn_Callback(hObject, eventdata, handles)
list = get(handles.fileListbox,'String');
file = list{get(handles.fileListbox,'Value')};
path = get(handles.pathBox,'String');

freqStart = get(handles.freqStartSlider, 'Value');
freqStop =  get(handles.freqStopSlider, 'Value');
crossFreq =  get(handles.crossFreqSlider, 'Value');
crossLength = get(handles.crossLengthSlider, 'Value');
nPoles1 = get(handles.nPoles1Slider, 'Value');
nPoles2 = get(handles.nPoles2Slider, 'Value');
lambda1 = get(handles.lambda1Slider, 'Value');
lambda2 = get(handles.lambda2Slider, 'Value');
ITER = str2num(get(handles.iterBox, 'String'));
Fs = str2double(handles.boxFs.String);
NFIR = str2double(handles.NFIRBox.String);
natSpline = handles.NAKSplineCheck.Value;

load([path '\' file]);  %fr,H
w = 2*pi*fr/Fs;
if handles.UseCppCheck.Value
    [Bm,Am,FIR] = filterDesignCpp(w,H, freqStart, freqStop, crossFreq, crossLength, nPoles1, nPoles2,lambda1, lambda2,NFIR,ITER,Fs,natSpline);
else
    [Bm,Am,FIR] = filterDesignMatlab(w,H, freqStart, freqStop, crossFreq, crossLength, nPoles1, nPoles2,lambda1, lambda2,NFIR,ITER,Fs,~natSpline);
end
    fname = [handles.exportCoeffsBox.String '.mat'];

if handles.exportParametersCheck.Value == 1
    save(fname,'Bm','Am','FIR','Fs','freqStart','freqStop','crossFreq','crossLength','nPoles1','nPoles2','lambda1','lambda2','ITER');
else
    if handles.CsvCheck.Value
        exportCsv(fname,Bm,Am,FIR);
    else
        save(fname,'Bm','Am','FIR');
    end
end


% --- Executes on button press in fixedNPolesCheck.
function fixedNPolesCheck_Callback(hObject, eventdata, handles)
if hObject.Value == 1
    handles.fixedNPolesBox.Enable = 'on';
    handles.nPoles2Slider.Enable = 'off';
    handles.nPoles2Box.Enable = 'off';
    fixedNPoles = str2double( handles.fixedNPolesBox.String);
    nPolesLow = handles.nPoles1Slider.Value;
    handles.nPoles2Box.String = num2str(fixedNPoles-nPolesLow);
    handles.nPoles2Slider.Value = fixedNPoles-nPolesLow;
    plotIt(handles);
else
    handles.fixedNPolesBox.Enable = 'off';
    handles.nPoles2Slider.Enable = 'on';
    handles.nPoles2Box.Enable = 'on';
end

function fixedNPolesBox_Callback(hObject, eventdata, handles)
fixedNPoles = str2double( hObject.String);
nPolesLow = handles.nPoles1Slider.Value;
set(handles.nPoles2Box, 'String', num2str(fixedNPoles-nPolesLow));
set(handles.nPoles2Slider, 'Value', fixedNPoles-nPolesLow);

plotIt(handles);

% --- Executes during object creation, after setting all properties.
function fixedNPolesBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in exportParametersCheck.
function exportParametersCheck_Callback(hObject, eventdata, handles)



function boxFs_Callback(hObject, eventdata, handles)
plotIt(handles);

% --- Executes during object creation, after setting all properties.
function boxFs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in viewOnlyCheck.
function viewOnlyCheck_Callback(hObject, eventdata, handles)
if hObject.Value == 1
    enab = 'off';
else   
    enab = 'on';
end
handles.crossFreqBox.Enable = enab;
handles.crossFreqSlider.Enable = enab;
handles.crossLengthBox.Enable = enab;
handles.crossLengthSlider.Enable = enab;
handles.nPoles1Box.Enable = enab;
handles.nPoles1Slider.Enable = enab;
handles.nPoles2Box.Enable = enab;
handles.nPoles2Slider.Enable = enab;
handles.lambda1Box.Enable = enab;
handles.lambda1Slider.Enable = enab;
handles.lambda2Box.Enable = enab;
handles.lambda2Slider.Enable = enab;
handles.NFIRBox.Enable = enab;
handles.iterBox.Enable = enab;
handles.exportCoeffsBtn.Enable = enab;
handles.fixedNPolesCheck.Enable = enab;


% --- Executes on button press in infoBtn.
function infoBtn_Callback(hObject, eventdata, handles)
Opt.Interpreter = 'tex';
Opt.WindowStyle = 'modal';
msgbox({'First enter path to location of stored files with frequency responses. Press \bfLoad\rm.';...
    'You should see available files in provided folder. Selecting different file will automatically create a new filter design and plot the result.';'';...
    '\bfFormat of files\rm - design mode';...
    'File should contain 2 fields:';'{\bffr} - column vector of frequencies in Hz';...
    '{\bfH} - column vector of frequency response (absolute, not in dB)'; '';...
    '\bfFormat of files\rm - view only mode';...
    'File should contain 3 fields: {\bfAm}, {\bfBm}, {\bfFIR} in the same format, that this app exports it.';...
    'It may containg fields target {\bffr} and {\bfH}, if so, response is also plotted. ';'';...
    '{\bfAlgorithm}';...
    'Filter design is based on Balazs Bank dual warping method. Frequency response is devided into 2 parts, each frequency vector is warped by lambda parameter and poles of filter are calculated using Kallman method.';...
    'Then iterative Steiglitz-McBride method optimizes the poles positions. Filter is divided into parallel second-order sections and zeros of filters are calculated.';...
    '';'You can see poles position on the bottom plot.';' ';...
    'Inside the algorithms spline interpolation is used. C++ uses linear extrapolation, MATLAB uses natural spline. Both can be switched to not-a-knot spline'; '';...
    'You can export the coefficients of parallel second-order sections and FIR part of the filter by {\bfexport} button. File will be located at current MATLAB directory.';...
    'Checking {\bfexport parameters} will export all parameters used in the design.';...
    '';'For more information see:';'http://home.mit.bme.hu/~bank/parfilt/';
    '';'App by Matou Vrbík';'matousvrbik@gmail.com'},...
    'How to use this app', 'help', Opt);

function NFIRBox_Callback(hObject, eventdata, handles)
plotIt(handles)

% --- Executes during object creation, after setting all properties.
function NFIRBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NAKSplineCheck.
function NAKSplineCheck_Callback(hObject, eventdata, handles)
plotIt(handles);


% --- Executes on button press in UseCppCheck.
function UseCppCheck_Callback(hObject, eventdata, handles)
plotIt(handles);


% --- Executes on button press in CsvCheck.
function CsvCheck_Callback(hObject, eventdata, handles)
if hObject.Value == 1
    handles.exportParametersCheck.Enable = 'off';
    handles.exportParametersCheck.Value = 0;
else
  handles.exportParametersCheck.Enable = 'on';
end
