function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 05-Jun-2017 19:48:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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

% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
addpath(fullfile(fileparts(mfilename('fullpath')),'Functions'));
addpath(fullfile(fileparts(mfilename('fullpath')),'Saved'));
set(handles.slider1, 'Value', 0);
set(handles.text19, 'String', '-');

% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mask_pp;
% Get default command line output from handles structure
varargout{1} = mask_pp;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes4,'reset');cla(handles.axes5,'reset');cla(handles.axes2,'reset');
cla(handles.axes6,'reset');cla(handles.axes7,'reset');cla(handles.axes9,'reset');
cla(handles.axes8,'reset');cla(handles.axes10,'reset');cla(handles.axes1,'reset');

set(handles.text18,'String','null');set(handles.text18,'Enable','off');
set(handles.text7,'String','0.0');set(handles.text7,'Enable','off');
global mask_pp;
global I;
[fn pn] = uigetfile({'*.jpg;*.png'},'Select Image File');

%nulis di kotak dialog lokasi nya
complete = strcat('Location : ',pn,fn); 
set(handles.edit1,'string',complete);

%buat baca gambar
I=imread([pn,fn]);

axes(handles.axes1);imshow(I,[]);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Eksekusi
set(handles.slider1, 'Value', 1);
global I;
global debug;
global mask_pp;
I_ori = I;
I=double(I);
set(handles.text19, 'String', '1');

% Simpen waktu
fprintf('Start ');
main_time = tic();
fprintf(' \n >')

%% Setting Parameter

lbr = size(I_ori,1);
tgg = size(I_ori,2);

param.PR.th_mask = 20;

% Parameters ZM
param.ZM.SZ         = 16;
param.ZM.ORDER      = 6;
param.ZM.radiusNum  = 26; 
param.ZM.anglesNum  = 32;

% Parameters Matching
param.M.minOFF          = 0.85;
param.M.minNNF          = 10; % jarak minimum kordi NNF
param.M.sizePatch       = 7;
param.M.midPatch        = floor(param.M.sizePatch/2)+1; % Nilai tengah patch
param.M.minItrNNF       = 0.9; %utk iter auto
param.M.iter_ini        = 10;

% Parameter Iterasi Auto
param.M.difer_mask      = 30; % perbedaaan MIN mask tiap iter, kalo dibawah ini iter berhenti
param.M.iter_auto_max   = 20;
param.M.iter_auto_min   = 2;

% Parameters Post Processing
param.PP.diameter       = 16;
param.PP.th2_dist2      = 50*50; % T^2_{D2} = minimum diatance between clones
param.PP.th2_dlf        = 300;   % T^2_{\epsion} = threshold on DLF error
param.PP.th_sizeA       = 5;  % T_{S} = minimum size of clones
param.PP.th_sizeB_pr    = 0.05;  
param.PP.th_sizeB       = round((lbr*tgg)*param.PP.th_sizeB_pr);
param.PP.rd_median      = 4;     % \rho_M = radius of median filter
param.PP.rd_dlf         = 6;     % \rho_N = radius of DLF patch
param.PP.rd_dil         = param.PP.rd_dlf + param.PP.rd_median; % \rho_D = radius for dilatetion

% Setting Iterasi
rd2=get(handles.radiobutton2,'Value');
if(rd2==1)
    %% iterasi patchmatch manual
    box_itr = get(handles.edit2,'String');
    itr = str2num(box_itr);
    param.M.main_iteration=itr;
    set(handles.text16,'String',num2str(itr));drawnow 
    
    %% Setting Slider
    set(handles.slider1, 'Min', 1);
    set(handles.slider1, 'Max', param.M.main_iteration);
    set(handles.slider1, 'Value', 1);drawnow 
    
    %% Ekstraksi Ciri: Zernike Moment
%     feat=I;
    feat = ZernikeMoment(I,param.ZM);

    %% Matching: PatchMatch
    feat  = (feat-min(feat(:)))./(max(feat(:))-min(feat(:))); %mPM requires the features to be in [0,1]
    [ debug, mask_pp ] = PM_DLF_GUI(feat(:,:,:,1), param,handles,I);
    
else 
    %% Iterasi Patchmatch Auto
    param.M.main_iteration=param.M.iter_auto_max;
    set(handles.slider1, 'Min', 1);
    set(handles.slider1, 'Max', param.M.main_iteration);drawnow 
    set(handles.slider1, 'Value', 1);
    set(handles.text16,'String',strcat('Max iterasi :',' ',num2str(param.M.main_iteration)));drawnow 
    
    %% Ekstraksi Ciri: Zernike Moment
    feat = ZernikeMoment(I,param.ZM);

    %% Matching: PatchMatch
    feat  = (feat-min(feat(:)))./(max(feat(:))-min(feat(:))); %mPM requires the features to be in [0,1]
    [ debug, mask_pp ] = PM_DLF_GUI_OTO(feat(:,:,:,1), param,handles,I);
end

fprintf('End ');
fprintf(' \n >')
save('debug_gui.mat','debug');
save('mask_pp_gui.mat','mask_pp');
main_time = toc(main_time)
set(handles.text7,'Enable','on');
set(handles.text7,'String',num2str(main_time));%time
% set(handles.text11,'String',num2str(main_time));%fmeasure
% set(handles.text14,'String',num2str(main_time));%Accuracy

varargout=debug;



% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton1.
function radiobutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function radiobutton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uipanel7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
manual = get(handles.radiobutton2,'Value');
if(manual == 1)
    set(handles.text5,'Enable','on');
    set(handles.edit2,'Enable','on');
else
    set(handles.text5,'Enable','off');
    set(handles.edit2,'Enable','off');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
%% slider
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% 
% global debug;
% global mask_pp;
load debug_gui;
load mask_pp_gui;
global I;
slide = get(handles.slider1, 'Value');drawnow

if (isreal(slide))
    slide=round(slide);
    set(handles.slider1,'Value',slide);drawnow
end
set(handles.text19, 'String', slide);

axes(handles.axes2);imshow(uint8(debug.NNF_prp_rs(slide).nnf(:,:,1)));
axes(handles.axes4);imshow(uint8(debug.NNF_prp_rs(slide).nnf(:,:,2)));drawnow

dist2 = MPFspacedist2((debug.NNF_prp_rs(slide).nnf(:,:,2)),debug.NNF_prp_rs(slide).nnf(:,:,1));
max_dist = max(sqrt(dist2(:)));
axes(handles.axes5);imshow(sqrt(dist2),[0, max_dist]); colormap(jet());

dist2 = MPFspacedist2(mask_pp.mpfy(slide).row,mask_pp.mpfx(slide).col);
axes(handles.axes6);imshow(sqrt(dist2),[0, max_dist]); colormap(jet()); %colorbar();drawnow 

DLF_db = 10*log10(mask_pp.pp_dlf(slide).pp_dlf.DLFerr); DLF_db(DLF_db<-50) = -50;
axes(handles.axes7);imshow(DLF_db,[]); colormap(jet());

axes(handles.axes9);imshow(double(repmat(mask_pp.pp_dlf(slide).pp_dlf.maskDLF,[1,1,3])));
axes(handles.axes8);imshow(double(repmat(mask_pp.pp_dlf(slide).pp_dlf.maskMPF,[1,1,3])));
axes(handles.axes10);imshow(double(repmat(mask_pp.mask(slide).pp_dlf,[1,1,3])));drawnow

mat_off=debug.mat_off_all(slide).disp_ofs;
mask=debug.mask(slide).mask;
mask2=debug.mask2(slide).mask2;
param_off=debug.param_off(slide).param_off;

axes(handles.axes1);imshow(uint8(I),[]);
if(max(max(mask2))~=0)
    hold on
    for i =param_off
       plot(mat_off(i,2), mat_off(i,1), 'ro');
       plot(mat_off(i,4), mat_off(i,3), 'ro');

       Y=[mat_off(i,1) mat_off(i,3)];
       X=[mat_off(i,2) mat_off(i,4)];
       line(X,Y,'Color','y');
    end
end
drawnow


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
