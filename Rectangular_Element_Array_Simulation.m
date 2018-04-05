function varargout = Rectangular_Element_Array_Simulation(varargin)
%% Rectangular_Element_Array_Simulation.m
%This script will execute a simulation of a phased transducer array for use in
%acoustoelectric cardiac imaging (ACI)
%Modified from Focus_Array_Pulse2.m written by Yexian Qin/Andy Tseng
%Written by Alexander Alvarez 10/27/17; v2.0 written on 10/30/17;
%clear;
%close all;
%clc;

%Set
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OpenFnc_Callback, ...
                   'gui_OutputFcn',  @OutputFnc_Callback, ...
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

function OpenFnc_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Construct the components of a new GUI
global TransParams
handles.output = hObject;
guidata(hObject, handles);
handles.p.Shape = 'Curved Strip';
handles.p.ElemX = 18;
handles.p.ElemY = 7;
handles.p.KerfX = 1e-4;
handles.p.KerfY = 1e-4;
handles.p.Width = 1.5e-3;
handles.p.Height = 13.3e-3;
handles.p.FocX = 0;
handles.p.FocY = 0;
handles.p.FocZ = 50e-3;
handles.p.Rad = 50e-3;
handles.p.Pressure = 0;
handles.p.CenFreq = 1.8e6;
handles.p.Power = 40;
handles.p.Efficiency = 0.75;
TransParams = handles.p;
reset(handles);

function reset(handles)
global TransParams
set(handles.Shape, 'String', {'Curved Strip', 'Planar', 'Cylindrical Section','Enclosed Cylindrical Section'});
set(handles.ElemX, 'String', num2str(TransParams.ElemX));
set(handles.ElemY, 'String', num2str(TransParams.ElemY));
set(handles.KerfX, 'String', num2str(TransParams.KerfX*1e3));
set(handles.KerfY, 'String', num2str(TransParams.KerfY*1e3));
set(handles.Width, 'String', num2str(TransParams.Width*1e3));
set(handles.Height, 'String', num2str(TransParams.Height*1e3));
set(handles.FocX, 'String', num2str(TransParams.FocX*1e3));
set(handles.FocY, 'String', num2str(TransParams.FocY*1e3));
set(handles.FocZ, 'String', num2str(TransParams.FocZ*1e3));
set(handles.Rad, 'String', num2str(TransParams.Rad*1e3));
set(handles.Pressure, 'String', num2str(TransParams.Pressure/1e6));
set(handles.CenFreq, 'String', num2str(TransParams.CenFreq/1e6));
set(handles.Power, 'String', num2str(TransParams.Power));
set(handles.Efficiency, 'String', num2str(TransParams.Efficiency));
% --- Executes on selection change in Shape.
function ArrayType_Callback(hObject, ~, ~)
% hObject    handle to Shape (see GCBO)
global TransParams
Array_Conts = cellstr(get(hObject,'String'));
TransParams.Shape = Array_Conts{get(hObject,'Value')};

% --- Executes during object creation, after setting all properties.
function ArrayType_CreateFcn(hObject, ~, ~)
% hObject    handle to Shape (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function XElem_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to XElem_Callback (see GCBO)
TransParams.ElemX = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function XElem_CreateFcn(hObject, ~, ~)
% hObject    handle to XElem_CreateFcn (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function YElem_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to ElemY (see GCBO)
TransParams.ElemY = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function YElem_CreateFcn(hObject, ~, ~)
% hObject    handle to ElemY (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function KerfX_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to KerfX (see GCBO)
TransParams.KerfX = (str2double(get(hObject,'String')))/1e3;

% --- Executes during object creation, after setting all properties.
function KerfX_CreateFcn(hObject, ~, ~)
% hObject    handle to KerfX (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function KerfY_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to KerfY (see GCBO)
TransParams.KerfY = (str2double(get(hObject,'String')))/1e3;

% --- Executes during object creation, after setting all properties.
function KerfY_CreateFcn(hObject, ~, ~)
% hObject    handle to KerfY (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Width_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to KerfX (see GCBO)
TransParams.Width = (str2double(get(hObject,'String')))/1e3;

% --- Executes during object creation, after setting all properties.
function Width_CreateFcn(hObject, ~, ~)
% hObject    handle to KerfX (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Height_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to KerfY (see GCBO)
TransParams.Height = (str2double(get(hObject,'String')))/1e3;

% --- Executes during object creation, after setting all properties.
function Height_CreateFcn(hObject, ~, ~)
% hObject    handle to KerfY (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RadCurv_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to KerfY (see GCBO)
TransParams.Rad = (str2double(get(hObject,'String')))/1e3;

% --- Executes during object creation, after setting all properties.
function RadCurv_CreateFcn(hObject, ~, ~)
% hObject    handle to KerfY (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FocX_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to KerfY (see GCBO)
TransParams.FocX = (str2double(get(hObject,'String')))/1e3;

% --- Executes during object creation, after setting all properties.
function FocX_CreateFcn(hObject, ~, ~)
% hObject    handle to KerfY (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function FocY_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to KerfY (see GCBO)
TransParams.FocY = (str2double(get(hObject,'String')))/1e3;

% --- Executes during object creation, after setting all properties.
function FocY_CreateFcn(hObject, ~, ~)
% hObject    handle to KerfY (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function FocZ_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to KerfY (see GCBO)
TransParams.FocZ = (str2double(get(hObject,'String')))/1e3;

% --- Executes during object creation, after setting all properties.
function FocZ_CreateFcn(hObject, ~, ~)
% hObject    handle to KerfY (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function CenFreq_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to KerfY (see GCBO)
TransParams.CenFreq = (str2double(get(hObject,'String')))*1e6;

% --- Executes during object creation, after setting all properties.
function CenFreq_CreateFcn(hObject, ~, ~)
% hObject    handle to KerfY (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Pressure_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to KerfY (see GCBO)
TransParams.Pressure = (str2double(get(hObject,'String')))*1e6;

% --- Executes during object creation, after setting all properties.
function Pressure_CreateFcn(hObject, ~, ~)
% hObject    handle to KerfY (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function InPower_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to KerfY (see GCBO)
TransParams.Power = (str2double(get(hObject,'String')));

% --- Executes during object creation, after setting all properties.
function InPower_CreateFcn(hObject, ~, ~)
% hObject    handle to KerfY (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Efficiency_Callback(hObject, ~, ~)
global TransParams
% hObject    handle to Efficiency (see GCBO)
TransParams.Efficiency = (str2double(get(hObject,'String')))*1e6;

% --- Executes during object creation, after setting all properties.
function Efficiency_CreateFcn(hObject, ~, ~)
% hObject    handle to Efficiency (see GCBO)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function varargout = OutputFnc_Callback(~, ~, handles)
varargout{1} = handles.output; 
% --- Executes on button press in PressureSimulate.

function PressureSimulate_Callback(~, ~, ~)
global TransParams
RecArray_SIM_forRuss_v2(TransParams);


% --- Executes on button press in AESimulate.
function AESimulate_Callback(~, ~, ~)
global TransParams
RecArray_SIM_v13(TransParams);
