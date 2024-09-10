function varargout = GravityTie(varargin)
% GRAVITYTIE MATLAB code for GravityTie.fig
%      GRAVITYTIE, by itself, creates a new GRAVITYTIE or raises the existing
%      singleton*.
%
%      H = GRAVITYTIE returns the handle to a new GRAVITYTIE or the handle to
%      the existing singleton*.
%
%      GRAVITYTIE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRAVITYTIE.M with the given input arguments.
%
%      GRAVITYTIE('Property','Value',...) creates a new GRAVITYTIE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GravityTie_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GravityTie_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GravityTie

% Last Modified by GUIDE v2.5 17-May-2021 08:09:42

% Mods by MQT v.3 15-April-2024 onboard SKQ202406T 
% nanmean --> replaced by mean (tested on Matlab 2024b windows)
% Sikuliaq (grav_dgs_33_proc....Z) and Thompson (ATM_GRAV_PROC_.RAW) have diff data format for File Data 
% Also made mods for better gui navigation - removing non-essential
% functions.

% Begin initialization code - DO NOT EDIT --------------------------------
% ------------------------------------------------------------------------
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GravityTie_OpeningFcn, ...
                   'gui_OutputFcn',  @GravityTie_OutputFcn, ...
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
% End initialization code - DO NOT EDIT-----------------------------------
%
% --- Executes just before GravityTie is made visible.
function GravityTie_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GravityTie (see VARARGIN)

% Choose default command line output for GravityTie
tf=ispc;
if isequal(tf,0)
    if ismac || isunix
    	tf='unix';
        handles.os='/';
    end
else
    tf='pc';
    handles.os='\';
end


handles.output = hObject;
set(findall(handles.uibuttongroup_A1, '-property', 'enable'), 'enable', 'off');
set(findall(handles.uibuttongroup_B, '-property', 'enable'), 'enable', 'off');
set(findall(handles.uibuttongroup_A2, '-property', 'enable'), 'enable', 'off');
set(findall(handles.report, '-property', 'enable'), 'enable', 'off');
set(findall(handles.panel_dgs, '-property', 'enable'), 'enable', 'off');
set(findall(handles.height, '-property', 'enable'), 'enable', 'off');

set(handles.uibuttongroup_A1,'ForegroundColor',[0.5 0.5 0.5]);
set(handles.uibuttongroup_A2,'ForegroundColor',[0.5 0.5 0.5]);
set(handles.uibuttongroup_B,'ForegroundColor',[0.5 0.5 0.5]);
set(handles.report,'ForegroundColor',[0.5 0.5 0.5]);
set(handles.panel_dgs,'ForegroundColor',[0.5 0.5 0.5]);
set(handles.height,'ForegroundColor',[0.5 0.5 0.5]);

set(handles.stationNo,'Enable','off');
set(handles.absgrav,'Enable','off');
set(handles.landmeter,'Enable','off');
set(handles.loadDgs,'Enable','off');
set(handles.plotdgs,'Enable','off')
set(handles.loadTie,'Enable','on')
set(handles.axes2,'Box','on');

lm_db_file = ['database' handles.os 'landmeters.db'];
% lm_db_file = 'landmeters.db';
lm = loadLandMeter(lm_db_file);
lmString{1} = 'Select a land meter';
lmString{2} = 'Other';
for i = 1:length(lm)
    lmString{i+2} = lm{i};
end
handles.nLandmeter = i;
set(handles.landmeter,'String',lmString);



ship_db_file = ['database' handles.os 'ships.db'];
% ship_db_file = 'ships.db';
ships = loadShip(ship_db_file);
shipString{1} = 'Select a ship';
for i = 1:length(ships)
    shipString{i+1} = ships{i};
end
set(handles.shipmenu,'String',shipString);

tie_db_file = ['database' handles.os 'stations.db'];
% tie_db_file = 'stations.db';
handles.allstations = loadGravitySta(tie_db_file);

stationString{1} = 'Select a station';
for i = 1:handles.allstations.nStation
    stationString{i+1} = handles.allstations.longName{i};
end

set(handles.stationmenu,'String',stationString);

    
handles.faafactor = 0.3086;
handles.d2s = 24*3600;
handles.ifloadTie = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GravityTie wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GravityTie_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
guidata(hObject, handles);


% --- Executes on selection change in shipmenu.
function shipmenu_Callback(hObject, eventdata, handles)
% hObject    handle to shipmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns shipmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from shipmenu
if handles.iflandtie.Value == 1 && handles.landmeter.Value == 1
    f = errordlg('Missing land meter #','Selection Error');
end
val = get(handles.shipmenu,'Value');
handles.ship = handles.shipmenu.String{val};
set(handles.loadTie,'Enable','off');
guidata(hObject, handles);



% --- Executes on selection change in stationmenu.
function stationmenu_Callback(hObject, eventdata, handles)
% hObject    handle to stationmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns stationmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stationmenu
val = get(handles.stationmenu,'Value');
handles.B.StaIndex = val;
handles.B.StaName = handles.stationmenu.String{val};
handles.B.StaNo = handles.allstations.Number{val-1};
handles.B.knownGrav = handles.allstations.Grav{val-1};

% set(handles.stationNo,'enable','on');
% set(handles.absgrav,'enable','on');
if contains(handles.B.StaName,'Other')
    set(handles.stationNo,'enable','on');
    set(handles.absgrav,'enable','on');
%     set(handles.stationNo,'String','');
%     set(handles.absgrav,'String','');
    
    

else
    set(handles.stationNo,'String',handles.B.StaNo);
    set(handles.absgrav,'String',handles.B.knownGrav);
end
% handles.B.StaNo = handles.stationNo.String;
% handles.B.knownGrav = handles.absgrav.String;


if handles.iflandtie.Value
    set(handles.uibuttongroup_A1,'ForegroundColor','k');
    set(handles.uibuttongroup_A2,'ForegroundColor','k');
    set(handles.uibuttongroup_B,'ForegroundColor','k');
    set(findall(handles.uibuttongroup_A1, '-property', 'enable'), 'enable', 'on');
    set(findall(handles.uibuttongroup_B, '-property', 'enable'), 'enable', 'on');
    set(findall(handles.uibuttongroup_A2, '-property', 'enable'), 'enable', 'on');
    
%     set(handles.computelandtie,'Enable','on');
    set(handles.newlandtietest,'Enable','off');
    set(handles.landtievalue,'Enable','off');
%     set(handles.A1Time1,'Enable','off');
%     set(handles.A1Time2,'Enable','off');
%     set(handles.A1Time3,'Enable','off');
%     set(handles.B1Time1,'Enable','off');
%     set(handles.B1Time2,'Enable','off');
%     set(handles.B1Time3,'Enable','off');
%     set(handles.A2Time1,'Enable','off');
%     set(handles.A2Time2,'Enable','off');
%     set(handles.A2Time3,'Enable','off');
end
set(handles.height,'ForegroundColor','k');
set(findall(handles.height, '-property', 'enable'), 'enable', 'on');
%     set(handles.heightTime1,'Enable','off');
%     set(handles.heightTime2,'Enable','off');
%     set(handles.heightTime3,'Enable','off');
% end

guidata(hObject, handles);


% --- Executes on button press in iflandtie.
function iflandtie_Callback(hObject, eventdata, handles)
% hObject    handle to iflandtie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of iflandtie
if handles.iflandtie.Value
    % land tie
    set(handles.landmeter,'Enable','on');
    set(handles.loadTie,'Enable','off');
    
else
end
guidata(hObject, handles);

function landmeter_Callback(hObject, eventdata, handles)
% hObject    handle to landmeter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of landmeter as text
%        str2double(get(hObject,'String')) returns contents of landmeter as a double
val = get(handles.landmeter,'Value');
if val == 2
    dlgTitle = '';
    dims = [1 100];
    prompt = {'\fontsize{12} Input land meter name: '};
    def = {''};
    AddOpts.Resize='on';
    AddOpts.Interpreter='tex';
    AddOpts.WindowStyle='normal';
    answer = inputdlg(prompt,dlgTitle,dims,def,AddOpts);
    SN = answer{1};
    
%     
%     % Add the new meter to the landmeters database
%     lm_db_file = ['database' handles.os 'landmeters.db'];
%     handles.nLandmeter = handles.nLandmeter+1;
%     [fid,msg] = fopen(lm_db_file,'a+');
%     if fid <0
%         error('Failed to open filename "%s" because: "%s"\n',lm_db_file,msg);
%     else
%         fprintf(fid,['[LAND_METER_' num2str(handles.nLandmeter) ']\n']);
%         fprintf(fid,['SN=' SN]);
%         fprintf(fid,['\nTABLE=' SN '.CAL']);
%         fprintf(fid,'\n\n');
%     end
%     fclose(fid);
end
if val > 2
    SN = handles.landmeter.String{val};
end
handles.calFile = ['database' handles.os 'land-cal' handles.os SN '.CAL'];
% handles.calFile = [SN '.CAL'];

handles.landmeterSN = SN;
guidata(hObject, handles);

% --- Executes on button press in computelandtie.
function computelandtie_Callback(hObject, eventdata, handles)
% hObject    handle to computelandtie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.B.StaNo = handles.stationNo.String;
handles.B.knownGrav = handles.absgrav.String;
for i = 1:3
    A1Readings(i) = str2double(handles.(['A1Read' num2str(i)]).String);
    A2Readings(i) = str2double(handles.(['A2Read' num2str(i)]).String);
    B1Readings(i) = str2double(handles.(['B1Read' num2str(i)]).String);
    
    A1Time(i) = handles.d2s*datenum(handles.(['A1Time' num2str(i)]).String);
    A2Time(i) = handles.d2s*datenum(handles.(['A2Time' num2str(i)]).String);
    B1Time(i) = handles.d2s*datenum(handles.(['B1Time' num2str(i)]).String);
    
end
if isfield(handles,'calTable')
    calTable = handles.calTable;
elseif isfile(handles.calFile)
    fid = fopen(handles.calFile,'r');
    Cols = textscan(fid,'%f%f%f','HeaderLines',1);
    fclose(fid);
    for i = 1:length(Cols{1,1})
        calTable(i,:) = [Cols{1,1}(i) Cols{1,2}(i) Cols{1,3}(i)];
    end
else
    calTable = [];
end


[mgalsA1,calTable,calTableupdate] = getLandTable(calTable,A1Readings);
[mgalsB1,calTable,calTableupdate] = getLandTable(calTable,B1Readings);
[mgalsA2,calTable,calTableupdate] = getLandTable(calTable,A2Readings);
handles.calTable = calTable;

% if calTableupdate == 1
% % Create the .CAL file for the new meter.
%     formatOut = 'yyyy/mm/dd';
%     t = datestr(now,formatOut)
%     fid = fopen(handles.calFile,'w+');
%     fprintf(fid,t);
%     fprintf(fid,'\n');
%     for i = 1:length(calTable(:,1))
%         fprintf(fid,'%04d\t%07.2f\t%7.5f\n',calTable(i,1),calTable(i,2),calTable(i,3));
%     end
%     
%     fclose(fid);
% end
   

avgTimeA1 = mean(A1Time);
avgTimeA2 = mean(A2Time);
avgTimeB1 = mean(B1Time);


avgMgalsA1 = mean(mgalsA1);
avgMgalsA2 = mean(mgalsA2);
avgMgalsB1 = mean(mgalsB1);

% Now do the calculation

aaTimeDelta = avgTimeA2 - avgTimeA1;
abTimeDelta = avgTimeB1 - avgTimeA1;
drift = (avgMgalsA2-avgMgalsA1)/aaTimeDelta;

dcAvgMgalsB1 = avgMgalsB1 - abTimeDelta * drift;
dcAvgMgalsA = avgMgalsA1;
% dcAvgMgalsA2 = avgMgalsA2 - aaTimeDelta * drift;
referenceG = str2double(handles.B.knownGrav);

gdiff = dcAvgMgalsA - dcAvgMgalsB1;
gravA = referenceG + gdiff;
handles.A.grav = gravA;
handles.A.drift = drift;
handles.A.avgMgals1 = avgMgalsA1;
handles.A.avgMgals2 = avgMgalsA2;
handles.B.avgMgals = avgMgalsB1;

handles.A.avgTime1 = avgTimeA1;
handles.A.avgTime2 = avgTimeA2;
handles.B.avgTime = avgTimeB1;

handles.A.dcAvgMgalsB1 = dcAvgMgalsB1;
handles.A.dcAvgMgalsA = dcAvgMgalsA;
handles.A.aaTimeDelta = aaTimeDelta;
handles.A.abTimeDelta = abTimeDelta;

set(handles.landtievalue,'String',num2str(gravA));
set(findall(handles.height, '-property', 'enable'), 'enable', 'on');
set(handles.height,'ForegroundColor',[0 0 0]);

set(handles.heightTime1,'Enable','on');
set(handles.heightTime2,'Enable','on');
set(handles.heightTime3,'Enable','on');
set(handles.newlandtietest,'Enable','on');
guidata(hObject, handles);




% --- Executes on button press in newlandtietest.
function newlandtietest_Callback(hObject, eventdata, handles)
% hObject    handle to newlandtietest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(handles.newlandtietest,'Value');
if val == 1
    
    set(handles.A1Read1,'String',' ');
    set(handles.A1Read2,'String',' ');
    set(handles.A1Read3,'String',' ');
    set(handles.A2Read1,'String',' ');
    set(handles.A2Read2,'String',' ');
    set(handles.A2Read3,'String',' ');
    set(handles.B1Read1,'String',' ');
    set(handles.B1Read2,'String',' ');
    set(handles.B1Read3,'String',' ');    
    set(handles.landtievalue,'String',' ');
    
    set(handles.A1Time1,'String',' ');
    set(handles.A1Time2,'String',' ');
    set(handles.A1Time3,'String',' ');
    set(handles.B1Time1,'String',' ');
    set(handles.B1Time2,'String',' ');
    set(handles.B1Time3,'String',' ');
    set(handles.A2Time1,'String',' ');
    set(handles.A2Time2,'String',' ');
    set(handles.A2Time3,'String',' ');


%     computelandtie_Callback(@computelandtie_Callback,eventdata, handles);
end
guidata(hObject, handles);

% --- Executes on button press in loadDgs.
function loadDgs_Callback(hObject, eventdata, handles)
% hObject    handle to loadDgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes2);
set(handles.axes2,'Visible','on');
axes(handles.axes2);
[fname, fpath] = uigetfile('*','MultiSelect','on');

ship = handles.ship;
if contains(lower(fname),'.dat') || contains(lower(fname),'proc')
    dat = loadDatDGS(ship,fname,fpath);    
else
    dat = loadRawDGS(ship,fname,fpath);
end

ts1 = datetime(handles.heightTime1.String,'InputFormat','yyyy/MM/dd HH:mm:ss');
t1 = convertTo(ts1,'posixtime');
ts2 = datetime(handles.heightTime3.String,'InputFormat','yyyy/MM/dd HH:mm:ss');
t2 = convertTo(ts2,'posixtime');

% dat.secTime = datenum(dat.yr,dat.mon,day,hh,mm,ss)*handles.d2s;

if min(dat.secTime) >= t2 | max(dat.secTime) <= t1
    f = errordlg('Data file doen''t match the water height measuring time. Select another file.',...
        'File Error');
else
    [~,h1] = min(abs(dat.secTime - t1));
    [~,h2] = min(abs(dat.secTime - t2));
    handles.meter.time = dat.secTime(h1:h2);
    handles.meter.rgrav = dat.rgrav(h1:h2);  
%     handles.meter.qcgrav = dat.qcgrav(h1:h2);

    set(handles.plotdgs,'Enable','on');
end
guidata(hObject, handles);

% --- Executes on button press in plotdgs.
function plotdgs_Callback(hObject, eventdata, handles)
% hObject    handle to plotdgs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
time = handles.meter.time;
dat = handles.meter.rgrav;
% dat = handles.meter.qcgrav;

% n = length(dat);
cla(handles.axes2);
set(handles.axes2,'Visible','on');
axes(handles.axes2);

ndata = length(time);
% filtertime = round(ndata/10);
filtertime = 361;
% sampling = 1;
% Taps=2*filtertime; 
% freq=1/filtertime;
% niquistf=sampling/2;
% wn=freq/niquistf;
% B = fir1(Taps,wn,blackman(Taps+1));
% fdat=filtfilt(B,1,dat);
[~,fdat] = gaussfilt(time,dat,filtertime);

% [tfilt,xfilt] = gaussfilt(time,dat,filterlength);
plot(dat(filtertime+1:end),'k','linewidth',1);
hold on
plot(fdat(filtertime+1:end),'r','linewidth',1)
% plot(tfilt(filterlength+1:end),xfilt(filterlength+1:end),'r','linewidth',1)
xlabel('Points');
ylabel('Grav (mGal)');
legend('Raw','Filtered')

handles.meter.avgGrav = mean(fdat(filtertime+1:end));

set(handles.computebias,'Enable','on')
% set(handles.loadDgs,'Enable','off')

guidata(hObject, handles);


% --- Executes on button press in computebias.
function computebias_Callback(hObject, eventdata, handles)
% hObject    handle to computebias (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set(handles.plotdgs,'Enable','off')
for i = 1:3
    heightPier(i)=str2double(handles.(['heightPier' num2str(i)]).String);

end
avgHeightPier = mean(heightPier);
handles.avgHeightPier = avgHeightPier;

set(findall(handles.report, '-property', 'enable'), 'enable', 'on');
set(handles.report,'ForegroundColor',[0,0,0]);

val = handles.iflandtie.Value;
if val
    pierGrav = str2double(handles.landtievalue.String);
    set(handles.heightA,'String',num2str(abs(avgHeightPier)));
else
    pierGrav = str2double(handles.absgrav.String);
    set(handles.meterTemp,'Enable','off');
    set(handles.latAdeg,'Enable','off');
    set(handles.latAmin,'Enable','off');
    set(handles.latAsec,'Enable','off');
    set(handles.lonAdeg,'Enable','off');
    set(handles.lonAmin,'Enable','off');
    set(handles.lonAsec,'Enable','off');
    set(handles.stationAname,'Enable','off');
    set(handles.heightA,'Enable','off');
end
waterGrav = pierGrav + handles.faafactor*avgHeightPier;
bias = waterGrav - handles.meter.avgGrav;

handles.meter.bias = bias;
handles.waterGrav = waterGrav;

set(handles.biasvalue,'String',num2str(bias));

guidata(hObject, handles);




% --- Executes on button press in makeReport.
function makeReport_Callback(hObject, eventdata, handles)
% hObject    handle to makeReport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
     
% 
% t = datetime('now','Timezone', 'UTC');
% ts = datestr(t,'yyyy/mm/dd HH:MM:SS');

if contains(handles.B.StaName,'Other')
    handles.B.StaNo = handles.stationNo.String;
    handles.B.knownGrav = handles.absgrav.String;
end

[outFile,outPath] = uiputfile('GravityReport.txt'); 

fid = fopen([outPath handles.os outFile],'w+');
fprintf(fid,['\nShip: ', handles.ship]);
fprintf(fid,['\nPersonnel: ', handles.personnel.String]);
fprintf(fid,'\n\n#Base station:');
fprintf(fid,['\nName: ', handles.B.StaName]);
fprintf(fid,['\nNumber: ', handles.B.StaNo]);
fprintf(fid,['\nKnown absolute gravity (mGal): ', handles.B.knownGrav]);
fprintf(fid,'\n\n------------------- Land tie ----------------');
val = handles.iflandtie.Value;
if val
    if handles.ifloadTie == 0
        lval = get(handles.landmeter,'Value');
        SN = handles.landmeter.String{lval};
    else
        SN = handles.landmeter.String;
    end
    fprintf(fid,['\n\nLand meter #: ' SN]);
    fprintf(fid,['\nMeter temperature (deg): ' handles.meterTemp.String]);
    fprintf(fid,'\n\n#New station A: ');
    latA = dms2dd(handles.latAdeg.String,handles.latAmin.String,handles.latAsec.String);
    lonA = dms2dd(handles.lonAdeg.String,handles.lonAmin.String,handles.lonAsec.String);
    fprintf(fid,['\nName: ', handles.stationAname.String]);
    fprintf(fid,['\nLatitude (deg): ', num2str(latA)]);
    fprintf(fid,['\nLongitude (deg): ', num2str(lonA)]);
    fprintf(fid,['\nElevation (m): ', handles.heightA.String]);
    
    ts = datenum2datetime(handles.A.avgTime1);
    fprintf(fid,['\nUTC time and meter gravity (mGal) at A1: ' ts '  ' num2str(handles.A.avgMgals1)]);
    ts = datenum2datetime(handles.B.avgTime);
    fprintf(fid,['\nUTC time and meter gravity (mGal) at B: ', ts '  ' num2str(handles.B.avgMgals)]);
    ts = datenum2datetime(handles.A.avgTime2);
    fprintf(fid,['\nUTC time and meter gravity (mGal) at A2: ', ts '  ' num2str(handles.A.avgMgals2)]);
    
    avgMgalsA2 = num2str(handles.A.avgMgals2);
    avgMgalsA1 = num2str(handles.A.avgMgals1);
    avgMgalsB = num2str(handles.B.avgMgals);
    dtaa = num2str(handles.A.aaTimeDelta);
    dtab = num2str(handles.A.abTimeDelta);
    drift = num2str(handles.A.drift);
    dcAvgMgalsB = num2str(handles.A.dcAvgMgalsB1);
    dcAvgMgalsA = num2str(handles.A.dcAvgMgalsA);
    referenceG = handles.B.knownGrav;
    gravPier = num2str(handles.A.grav);
    
    fprintf(fid,['\nDelta_T_ab (s): ' dtab]);
    fprintf(fid,['\nDelta_T_aa (s): ' dtaa]);
    fprintf(fid,['\nDrift (mGal): (', avgMgalsA2 ' - ' avgMgalsA1 ')/' dtaa ' = ' drift]);
    fprintf(fid,['\nDrift corrected meter gravity at B (mGal): ', ...
        avgMgalsB ' - ' dtab ' * ' drift ' = ' dcAvgMgalsB]);
    fprintf(fid,['\nGravity at pier (mGal): ', ...
        referenceG ' + ' dcAvgMgalsA ' - ' dcAvgMgalsB ' = ' gravPier]);
    
    fprintf(fid,'\n\n------------------- End of land tie ----------------');
    
else
    fprintf(fid,'\n\nN/A');
    fprintf(fid,'\n\n------------------- End of land tie ----------------');
    gravPier = handles.B.knownGrav;
    fprintf(fid,['\n\nGravity at pier (mGal): ' gravPier]);
end
t = datenum(handles.heightTime1.String)*handles.d2s;
ts = datenum2datetime(t);
fprintf(fid,['\nUTC time and water height to pier (m) 1: ' ts '  ' handles.heightPier1.String]);
t = datenum(handles.heightTime2.String)*handles.d2s;
ts = datenum2datetime(t);
fprintf(fid,['\nUTC time and water height to pier (m) 2: ', ts '  ' handles.heightPier2.String]);
t = datenum(handles.heightTime3.String)*handles.d2s;
ts = datenum2datetime(t);
fprintf(fid,['\nUTC time and water height to pier (m) 3: ', ts '  ' handles.heightPier3.String]);
   


avgHeightPier = num2str(handles.avgHeightPier);
waterGrav = num2str(handles.waterGrav);
meterAvgGrav = num2str(handles.meter.avgGrav);
meterBias = num2str(handles.meter.bias);
    
fprintf(fid,['\n\nDgS meter gravity (mGal): ', meterAvgGrav]);
fprintf(fid,['\nAverage water height to pier (m): ', avgHeightPier]);
fprintf(fid,['\nGravity at water line (mGal): ' gravPier ...
    ' + 0.3086 * ' avgHeightPier ' = ' waterGrav]);
fprintf(fid,['\nDgS meter bias (mGal): ' waterGrav ...
    ' - ' meterAvgGrav ' = ' meterBias]);


fprintf(fid,'\n\n------------------- End ----------------');
fclose(fid);

f = msgbox('Report generated.');


function personnel_Callback(hObject, eventdata, handles)
% hObject    handle to personnel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of personnel as text
%        str2double(get(hObject,'String')) returns contents of personnel as a double


% --- Executes during object creation, after setting all properties.
function personnel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to personnel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function heightMeter1_Callback(hObject, eventdata, handles)
% hObject    handle to heightMeter1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of heightMeter1 as text
%        str2double(get(hObject,'String')) returns contents of heightMeter1 as a double


% --- Executes during object creation, after setting all properties.
function heightMeter1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heightMeter1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function heightMeter2_Callback(hObject, eventdata, handles)
% hObject    handle to heightMeter2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of heightMeter2 as text
%        str2double(get(hObject,'String')) returns contents of heightMeter2 as a double


% --- Executes during object creation, after setting all properties.
function heightMeter2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heightMeter2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function heightMeter3_Callback(hObject, eventdata, handles)
% hObject    handle to heightMeter3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of heightMeter3 as text
%        str2double(get(hObject,'String')) returns contents of heightMeter3 as a double


% --- Executes during object creation, after setting all properties.
function heightMeter3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heightMeter3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function heightA_Callback(hObject, eventdata, handles)
% hObject    handle to heightA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of heightA as text
%        str2double(get(hObject,'String')) returns contents of heightA as a double


% --- Executes during object creation, after setting all properties.
function heightA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heightA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function meterTemp_Callback(hObject, eventdata, handles)
% hObject    handle to meterTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of meterTemp as text
%        str2double(get(hObject,'String')) returns contents of meterTemp as a double


% --- Executes during object creation, after setting all properties.
function meterTemp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meterTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function A2Time1_Callback(hObject, eventdata, handles)
% hObject    handle to A2Time1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A2Time1 as text
%        str2double(get(hObject,'String')) returns contents of A2Time1 as a double


% --- Executes during object creation, after setting all properties.
function A2Time1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A2Time1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A2Read1_Callback(hObject, eventdata, handles)
% hObject    handle to A2Read1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A2Read1 as text
%        str2double(get(hObject,'String')) returns contents of A2Read1 as a double
t = datetime('now','Timezone', 'UTC');
ts = datestr(t,'yyyy/mm/dd HH:MM:SS');
set(handles.A2Time1,'String',ts);



% --- Executes during object creation, after setting all properties.
function A2Read1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A2Read1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A2Time2_Callback(hObject, eventdata, handles)
% hObject    handle to A2Time2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A2Time2 as text
%        str2double(get(hObject,'String')) returns contents of A2Time2 as a double


% --- Executes during object creation, after setting all properties.
function A2Time2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A2Time2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A2Read2_Callback(hObject, eventdata, handles)
% hObject    handle to A2Read2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A2Read2 as text
%        str2double(get(hObject,'String')) returns contents of A2Read2 as a double
t = datetime('now','Timezone', 'UTC');
ts = datestr(t,'yyyy/mm/dd HH:MM:SS');
set(handles.A2Time2,'String',ts);


% --- Executes during object creation, after setting all properties.
function A2Read2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A2Read2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A2Time3_Callback(hObject, eventdata, handles)
% hObject    handle to A2Time3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A2Time3 as text
%        str2double(get(hObject,'String')) returns contents of A2Time3 as a double


% --- Executes during object creation, after setting all properties.
function A2Time3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A2Time3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A2Read3_Callback(hObject, eventdata, handles)
% hObject    handle to A2Read3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A2Read3 as text
%        str2double(get(hObject,'String')) returns contents of A2Read3 as a double
t = datetime('now','Timezone', 'UTC');
ts = datestr(t,'yyyy/mm/dd HH:MM:SS');
set(handles.A2Time3,'String',ts);


% --- Executes during object creation, after setting all properties.
function A2Read3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A2Read3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B1Read3_Callback(hObject, eventdata, handles)
% hObject    handle to B1Read3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B1Read3 as text
%        str2double(get(hObject,'String')) returns contents of B1Read3 as a double
t = datetime('now','Timezone', 'UTC');
ts = datestr(t,'yyyy/mm/dd HH:MM:SS');
set(handles.B1Time3,'String',ts);



% --- Executes during object creation, after setting all properties.
function B1Read3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B1Read3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B1Time3_Callback(hObject, eventdata, handles)
% hObject    handle to B1Time3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B1Time3 as text
%        str2double(get(hObject,'String')) returns contents of B1Time3 as a double


% --- Executes during object creation, after setting all properties.
function B1Time3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B1Time3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B1Read2_Callback(hObject, eventdata, handles)
% hObject    handle to B1Read2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B1Read2 as text
%        str2double(get(hObject,'String')) returns contents of B1Read2 as a double
t = datetime('now','Timezone', 'UTC');
ts = datestr(t,'yyyy/mm/dd HH:MM:SS');
set(handles.B1Time2,'String',ts);



% --- Executes during object creation, after setting all properties.
function B1Read2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B1Read2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B1Time2_Callback(hObject, eventdata, handles)
% hObject    handle to B1Time2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B1Time2 as text
%        str2double(get(hObject,'String')) returns contents of B1Time2 as a double


% --- Executes during object creation, after setting all properties.
function B1Time2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B1Time2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B1Read1_Callback(hObject, eventdata, handles)
% hObject    handle to B1Read1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B1Read1 as text
%        str2double(get(hObject,'String')) returns contents of B1Read1 as a double
t = datetime('now','Timezone', 'UTC');
ts = datestr(t,'yyyy/mm/dd HH:MM:SS');
set(handles.B1Time1,'String',ts);


% --- Executes during object creation, after setting all properties.
function B1Read1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B1Read1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B1Time1_Callback(hObject, eventdata, handles)
% hObject    handle to B1Time1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B1Time1 as text
%        str2double(get(hObject,'String')) returns contents of B1Time1 as a double


% --- Executes during object creation, after setting all properties.
function B1Time1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B1Time1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A1Time1_Callback(hObject, eventdata, handles)
% hObject    handle to A1Time1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A1Time1 as text
%        str2double(get(hObject,'String')) returns contents of A1Time1 as a double
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function A1Time1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A1Time1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A1Read1_Callback(hObject, eventdata, handles)
% hObject    handle to A1Read1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A1Read1 as text
%        str2double(get(hObject,'String')) returns contents of A1Read1 as a double

t = datetime('now','Timezone', 'UTC');
ts = datestr(t,'yyyy/mm/dd HH:MM:SS');
set(handles.A1Time1,'String',ts);
    
% --- Executes during object creation, after setting all properties.
function A1Read1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A1Read1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A1Time2_Callback(hObject, eventdata, handles)
% hObject    handle to A1Time2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A1Time2 as text
%        str2double(get(hObject,'String')) returns contents of A1Time2 as a double


% --- Executes during object creation, after setting all properties.
function A1Time2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A1Time2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A1Read2_Callback(hObject, eventdata, handles)
% hObject    handle to A1Read2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A1Read2 as text
%        str2double(get(hObject,'String')) returns contents of A1Read2 as a double

t = datetime('now','Timezone', 'UTC');
ts = datestr(t,'yyyy/mm/dd HH:MM:SS');
set(handles.A1Time2,'String',ts);


% --- Executes during object creation, after setting all properties.
function A1Read2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A1Read2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A1Time3_Callback(hObject, eventdata, handles)
% hObject    handle to A1Time3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A1Time3 as text
%        str2double(get(hObject,'String')) returns contents of A1Time3 as a double


% --- Executes during object creation, after setting all properties.
function A1Time3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A1Time3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A1Read3_Callback(hObject, eventdata, handles)
% hObject    handle to A1Read3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A1Read3 as text
%        str2double(get(hObject,'String')) returns contents of A1Read3 as a double
t = datetime('now','Timezone', 'UTC');
ts = datestr(t,'yyyy/mm/dd HH:MM:SS');
set(handles.A1Time3,'String',ts);


% --- Executes during object creation, after setting all properties.
function A1Read3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A1Read3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function stationmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stationmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function stationNo_Callback(hObject, eventdata, handles)
% hObject    handle to stationNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stationNo as text
%        str2double(get(hObject,'String')) returns contents of stationNo as a double


% --- Executes during object creation, after setting all properties.
function stationNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stationNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function absgrav_Callback(hObject, eventdata, handles)
% hObject    handle to absgrav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of absgrav as text
%        str2double(get(hObject,'String')) returns contents of absgrav as a double


% --- Executes during object creation, after setting all properties.
function absgrav_CreateFcn(hObject, eventdata, handles)
% hObject    handle to absgrav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function shipmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shipmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function heightPier1_Callback(hObject, eventdata, handles)
% hObject    handle to heightPier1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of heightPier1 as text
%        str2double(get(hObject,'String')) returns contents of heightPier1 as a double

t = datetime('now','Timezone', 'UTC');
ts = datestr(t,'yyyy/mm/dd HH:MM:SS');
set(handles.heightTime1,'String',ts);

% --- Executes during object creation, after setting all properties.
function heightPier1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heightPier1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function heightPier2_Callback(hObject, eventdata, handles)
% hObject    handle to heightPier2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of heightPier2 as text
%        str2double(get(hObject,'String')) returns contents of heightPier2 as a double
t = datetime('now','Timezone', 'UTC');
ts = datestr(t,'yyyy/mm/dd HH:MM:SS');
set(handles.heightTime2,'String',ts);


% --- Executes during object creation, after setting all properties.
function heightPier2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heightPier2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function heightPier3_Callback(hObject, eventdata, handles)
% hObject    handle to heightPier3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of heightPier3 as text
%        str2double(get(hObject,'String')) returns contents of heightPier3 as a double

t = datetime('now','Timezone', 'UTC');
ts = datestr(t,'yyyy/mm/dd HH:MM:SS');
set(handles.heightTime3,'String',ts);
set(findall(handles.panel_dgs, '-property', 'enable'), 'enable', 'on');
set(handles.panel_dgs,'ForegroundColor',[0 0 0]);
set(handles.plotdgs,'Enable','off');
set(handles.computebias,'Enable','off');
set(handles.biasvalue,'Enable','off');

% --- Executes during object creation, after setting all properties.
function heightPier3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heightPier3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function landmeter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to landmeter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function latAdeg_Callback(hObject, eventdata, handles)
% hObject    handle to latAdeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latAdeg as text
%        str2double(get(hObject,'String')) returns contents of latAdeg as a double


% --- Executes during object creation, after setting all properties.
function latAdeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latAdeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function latAmin_Callback(hObject, eventdata, handles)
% hObject    handle to latAmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latAmin as text
%        str2double(get(hObject,'String')) returns contents of latAmin as a double


% --- Executes during object creation, after setting all properties.
function latAmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latAmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function latAsec_Callback(hObject, eventdata, handles)
% hObject    handle to latAsec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latAsec as text
%        str2double(get(hObject,'String')) returns contents of latAsec as a double


% --- Executes during object creation, after setting all properties.
function latAsec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latAsec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lonAdeg_Callback(hObject, eventdata, handles)
% hObject    handle to lonAdeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lonAdeg as text
%        str2double(get(hObject,'String')) returns contents of lonAdeg as a double


% --- Executes during object creation, after setting all properties.
function lonAdeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lonAdeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lonAmin_Callback(hObject, eventdata, handles)
% hObject    handle to lonAmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lonAmin as text
%        str2double(get(hObject,'String')) returns contents of lonAmin as a double


% --- Executes during object creation, after setting all properties.
function lonAmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lonAmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lonAsec_Callback(hObject, eventdata, handles)
% hObject    handle to lonAsec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lonAsec as text
%        str2double(get(hObject,'String')) returns contents of lonAsec as a double


% --- Executes during object creation, after setting all properties.
function lonAsec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lonAsec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stationAname_Callback(hObject, eventdata, handles)
% hObject    handle to stationAname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stationAname as text
%        str2double(get(hObject,'String')) returns contents of stationAname as a double


% --- Executes during object creation, after setting all properties.
function stationAname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stationAname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function landtievalue_Callback(hObject, eventdata, handles)
% hObject    handle to landtievalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of landtievalue as text
%        str2double(get(hObject,'String')) returns contents of landtievalue as a double


% --- Executes during object creation, after setting all properties.
function landtievalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to landtievalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function biasvalue_Callback(hObject, eventdata, handles)
% hObject    handle to biasvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of biasvalue as text
%        str2double(get(hObject,'String')) returns contents of biasvalue as a double


% --- Executes during object creation, after setting all properties.
function biasvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to biasvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function heightTime2_Callback(hObject, eventdata, handles)
% hObject    handle to heightTime2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of heightTime2 as text
%        str2double(get(hObject,'String')) returns contents of heightTime2 as a double


% --- Executes during object creation, after setting all properties.
function heightTime2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heightTime2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function heightTime3_Callback(hObject, eventdata, handles)
% hObject    handle to heightTime3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of heightTime3 as text
%        str2double(get(hObject,'String')) returns contents of heightTime3 as a double


% --- Executes during object creation, after setting all properties.
function heightTime3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heightTime3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function heightTime1_Callback(hObject, eventdata, handles)
% hObject    handle to heightTime1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of heightTime1 as text
%        str2double(get(hObject,'String')) returns contents of heightTime1 as a double


% --- Executes during object creation, after setting all properties.
function heightTime1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heightTime1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in loadTie.
function loadTie_Callback(hObject, eventdata, handles)
% hObject    handle to loadTie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.stationmenu,'Enable','off');
set(handles.shipmenu,'Enable','off');
handles.ifloadTie = 1;
[fname, fpath] = uigetfile('*.mat');
load([fpath fname]);
handles.ship = tie.ship;
handles.B = tie.B;
handles.iflandtie.Value = tie.iflandtie.Value;
handles.landmeter.String = tie.landmeter;
handles.heightA.String = tie.heightA;
handles.shipmenu.String = handles.ship;
handles.stationmenu.String = handles.B.StaName;
handles.stationNo.String = handles.B.StaNo;
handles.absgrav.String = handles.B.knownGrav;
handles.landtievalue.String = tie.landtievalue;
handles.calTable = tie.calTable;
for i = 1:3
    handles.(['heightTime' num2str(i)]).String = tie.(['heightTime' num2str(i)]);
    handles.(['heightPier' num2str(i)]).String = tie.(['heightPier' num2str(i)]);
    
end
if sum(strcmp(fieldnames(tie), 'A')) == 1
    handles.A = tie.A;
    for i = 1:3
        handles.(['A1Time' num2str(i)]).String = tie.(['A1Time' num2str(i)]);
        handles.(['B1Time' num2str(i)]).String = tie.(['B1Time' num2str(i)]);
        handles.(['A2Time' num2str(i)]).String = tie.(['A2Time' num2str(i)]);
        handles.(['A1Read' num2str(i)]).String = tie.(['A1Read' num2str(i)]);
        handles.(['A2Read' num2str(i)]).String = tie.(['A2Read' num2str(i)]);
        handles.(['B1Read' num2str(i)]).String = tie.(['B1Read' num2str(i)]);
    end
end

set(handles.loadDgs,'Enable','on');
guidata(hObject, handles);

% --- Executes on button press in saveTie.
function saveTie_Callback(hObject, eventdata, handles)
% hObject    handle to loadTie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[outFile,outPath] = uiputfile('tieParameters.mat'); 
tie.ship = handles.ship;
tie.B= handles.B;
tie.iflandtie.Value = handles.iflandtie.Value;
% val = handles.landmeter.Value;
try
    tie.landmeter = handles.landmeterSN;
catch
    tie.landmeter = [];
end
tie.heightA = handles.heightA.String;
tie.landtievalue = handles.landtievalue.String;
try
    tie.calTable = handles.calTable;
catch 
    tie.calTable = [];
end

for i = 1:3
    tie.(['heightTime' num2str(i)]) = handles.(['heightTime' num2str(i)]).String;
    tie.(['heightPier' num2str(i)]) = handles.(['heightPier' num2str(i)]).String;
    
end
if sum(strcmp(fieldnames(handles), 'A')) == 1
    tie.A= handles.A;
    for i = 1:3
        tie.(['A1Time' num2str(i)]) = handles.(['A1Time' num2str(i)]).String;
        tie.(['B1Time' num2str(i)]) = handles.(['B1Time' num2str(i)]).String;
        tie.(['A2Time' num2str(i)]) = handles.(['A2Time' num2str(i)]).String;
        tie.(['A1Read' num2str(i)]) = handles.(['A1Read' num2str(i)]).String;
        tie.(['A2Read' num2str(i)]) = handles.(['A2Read' num2str(i)]).String;
        tie.(['B1Read' num2str(i)]) = handles.(['B1Read' num2str(i)]).String;
    end
end

save([outPath outFile],'tie');
set(handles.loadTie,'Enable','off');

f = msgbox('Tie parameters are saved.');


%----------------------------------------------------------------------
%         sub functions
%----------------------------------------------------------------------

function lm = loadLandMeter(lm_db_file)
% lm_db_file = 'database/landmeters.db';
fh = fopen(lm_db_file,'r');
N = 1;
while ~feof(fh)
    line = fgetl(fh);
    if (startsWith(line,'SN'))  % section found
        value = processline(line);                
        lm{N} = value;
        N = N+1;
    end        
end

% -------------------------------------
function ships = loadShip(ship_db_file)    
% ship_db_file = 'database/ships.db';
fh = fopen(ship_db_file,'r');
N = 1;
while ~feof(fh)
    line = fgetl(fh);
    if (startsWith(line,'SHIP'))  % section found
        value = processline(line);                
        ships{N} = strrep(value,'"','');
        N = N+1;
    end        
end

%------------------------------------
function S = loadGravitySta(tie_db_file)    
fh = fopen(tie_db_file,'r');
while ~feof(fh)
    line = fgetl(fh);
    if (startsWith(line,'[STATION'))  % section found
        N = str2double(line(10:end-1));
        for j = 1:12
            line = fgetl(fh);
            value{j} = processline(line);                
        end
        S.Country{N} = strrep(value{1},'"','');
        S.State{N} = strrep(value{2},'"','');
        S.City{N} = strrep(value{3},'"','');
        S.Name{N} = strrep(value{4},'"','');
        S.longName{N} = [S.Country{N} ' - ' S.State{N} ' - ' S.Name{N}];
        S.Number{N} = strrep(value{5},'"','');
        S.Date{N} = strrep(value{6},'"','');
        S.Lat{N} = strrep(value{7},'"','');
        S.Lon{N} = strrep(value{8},'"','');
        S.Grav{N} = strrep(value{9},'"','');
        S.Active{N} = strrep(value{10},'"','');
        S.Confidence{N} = strrep(value{11},'"','');
%         S.PDF{N} = strrep(value{12},'"','');
    end       
end
S.nStation = N;

%------------------------------------
function value = processline(line)
% Processes a line read from a db file and
A = strsplit(line,'=');
value = A{2};

%------------------------------------
function [g,calTable,calTableupdate] = getLandTable(calTable,reading)
% table = 'database/land-cal/G-70.CAL';
calTableupdate = 0;
if isempty(calTable)
    [calTable,ind] = inputCalTable(calTable,reading(1));
     calTableupdate = 1;
end
for i = 1:3
    if length(calTable(:,1)) == 1
        m = reading(i) - calTable(1,1);
        if m > 100 || m < 0
            [calTable,ind] = inputCalTable(calTable,reading(i));
            calTableupdate = 1;
        else
            ind = 1;
        end
    else
        [m,h] = min(abs(reading(i)-calTable(:,1)));
        if m < 100
            if reading(i)<calTable(h,1) 
                if reading(i) - calTable(h-1,1) <100
                    ind = h - 1;
                else
                    [calTable,ind] = inputCalTable(calTable,reading(i));
                    calTableupdate = 1;
                end
            else
                ind = h;
            end

        else % table is incomplete
            [calTable,ind] = inputCalTable(calTable,reading(i));
            calTableupdate = 1;
        end
    end
    readBracket = calTable(ind,1);
    factor = calTable(ind,3);
    offset = calTable(ind,2);
    residualReading = reading(i) - readBracket;
    residual = residualReading .* factor;
    g(i) = offset + residual;
end
    

% %-----------------------------------
function [currentCalTable,ind] = inputCalTable(currentCalTable,reading)
dlgTitle = 'Input land meter calibration table';
dims = [1 100];
prompt = {['\fontsize{12} Input Counter reading, Value in mGals, and Factor for interval (comma separated) for the reading of ' num2str(reading)]};
def = {''};
AddOpts.Resize='on';
AddOpts.Interpreter='tex';
AddOpts.WindowStyle='normal';
answer = inputdlg(prompt,dlgTitle,dims,def,AddOpts);
str = strsplit(answer{1},',');
calTableNew = str2double(str);
currentCalTable = [currentCalTable;calTableNew]; % append to currentCalTable
[currentCalTable,hh] = sortrows(currentCalTable,1);% sort table
[~,ind] = max(hh);


%-------------------------------------

function dd = dms2dd(degStr,minStr,secStr)
sg = sign(str2double(degStr));
d = str2double(degStr);
m = str2double(minStr);
s = str2double(secStr);    
dd = d + sg * (m/60 + s/3600);

% -------------------------------------
function ts = datenum2datetime(dnum)
t = datetime(dnum/24/3600,'ConvertFrom','datenum');
ts = datestr(t,'yyyy/mm/dd HH:MM:SS');

% % ---------------------------------------
    
function dat = loadDatDGS(ship,fname,fpath)
    if contains(ship,'Thompson')
        fmt = ['%2d/%2d/%4d%2d:%2d:%f' repmat('%f',1,19) '%*[^\n]'];
    elseif contains(ship,'Sikuliaq')
        fmt = ['%*s %f%*f%*f%*f%*f%*f%*f%*f%*f%f%f%f%f%f%f%*f%*f%*f%d%d%d%d%d%f%*f'];
    else
%         if contains(ship,'NBP') || contains(ship,'Atlantis') || contains(ship,'Revelle') || contains(ship,'Ride')
        % selected columns:
        % 1. uncorrected gravity
        % 2. ve,  3. vcc, 4. al, 5. ax
        % 6. lat, 7. lon
        % 8. year,  9. month,  10. day,  11. hour,  12. min,  13. sec
        fmt = '%*f%f%*f%*f%*f%*f%*f%*f%*f%*f%f%f%f%f%f%f%*f%*f%*f%d%d%d%d%d%f%*f'; 
    end
    %% Constant
    d2s = 24*3600; % day to secs
    fid=fopen(fullfile(fpath,fname),'r');
    if contains(ship,'Sikuliaq')
    opt = {'HeaderLines',27};
    C = textscan(fid,fmt,'delimiter',',',opt{:});
    else
    C = textscan(fid,fmt,'delimiter',',');
    end
    fclose(fid);
    [yr,mon,day,hh,mm,ss,rgrav,ve,vcc,al,ax,lat,lon] = parsingDatDGS(C,ship);
    dat.yr = yr;
    dat.mon = mon;
    dat.day = day;
    dat.hh = hh;
    dat.mm = mm;
    dat.ss = ss;
    dat.rgrav = rgrav;
    dat.vcc = -0.000029*vcc;
    dat.ve = 0.00001*ve;
    dat.al = 0.00001*al;
    dat.ax = 0.00001*ax;
    dat.lat = lat;
    dat.lon = lon;
    dat.ts = datetime(yr,mon,day,hh,mm,ss); % date string
    dat.secTime = convertTo(dat.ts,'posixtime'); % numerical date in seconds
    dat.numTime = dat.secTime/d2s; % numeric date in days
    

function [yr,mon,day,hh,mm,ss,rgrav,ve,vcc,al,ax,lat,lon] = parsingDatDGS(C,ship) 
    if contains(ship,'Thompson')
        yr = double(C{3});
        mon = double(C{1});
        day = double(C{2});
        hh = double(C{4});
        mm = double(C{5});
        ss = double(C{6});
        ss = round(ss);
        rgrav= double(C{8});
        ve = double(C{17});
        vcc = double(C{18});
        al = double(C{19});
        ax = double(C{20});
        lat = double(C{21});
        lon = double(C{22});
    else
%         if contains(ship,'NBP')|| contains(ship,'Atlantis')|| contains(ship,'Revelle')|| contains(ship,'Ride')
    
        % selected columns:
        % 1. uncorrected gravity
        % 2. ve,  3. vcc, 4. al, 5. ax
        % 6. lat, 7. lon
        % 8. year,  9. month,  10. day,  11. hour,  12. min,  13. sec
        rgrav= C{1};
        ve = C{2};
        vcc = C{3};
        al = C{4};
        ax = C{5};
        lat = C{6};
        lon = C{7};
        yr = double(C{8});
        mon = double(C{9});
        day = double(C{10});
        hh = double(C{11});
        mm = double(C{12});
        ss = double(C{13});
        ss = round(ss);
    end
    