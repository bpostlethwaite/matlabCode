function varargout = Wavelet_Denoise(varargin)
% WAVELET_DENOISE Application M-file for Wavelet_Denoise.fig
%    FIG = WAVELET_DENOISE launch Wavelet_Denoise GUI.
%    WAVELET_DENOISE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 11-Nov-2004 09:30:52

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function varargout = pushbutton4_Callback(h, eventdata, handles, varargin)

[filename, pathname] = uigetfile( ...
    {'*.wav', 'All wav-Files (*.wav)'; ...
        '*.*','All Files (*.*)'}, ...
    'Select Wave File');
% If "Cancel" is selected then return

File = fullfile(pathname,filename);
[x1 fs nbits]=wavread(File);
x2=x1(1:65536);
axes(handles.axes1)
plot(x2);
set(handles.axes1,'XMinorTick','on')
grid on
handles.X1=x2;
guidata(h, handles);
handles.X2=fs;
guidata(h, handles);

set(handles.edit1,'String',File)
uiwait(msgbox('File Loaded','File Loaded','modal'));


% --------------------------------------------------------------------
function varargout = Denoise_Callback(h, eventdata, handles, varargin)


f1 = get(handles.edit1,'String')
x=findstr(f1,'\')
if x==[]
[s,N] = makesig(f1);
n = randn(1,N);

f2 = str2double(get(handles.edit4,'String'));
if f2==2
h = 2; 
elseif f2==4
h = 4; 
elseif f2==6;
h = 6; 
elseif f2==8;
h = 8 ; 
elseif f2==10;
h = 10; 
else
h = 2; 
end  
x = s + n/h; 
f2 = str2double(get(handles.edit3,'String'))
if f2==4
h = daubcqf(4); 
elseif f2==6
h = daubcqf(6); 
elseif f2==8;
h = daubcqf(8); 
elseif f2==10;
h = daubcqf(10); 
else
h = daubcqf(4); 
end
f1 = str2double(get(handles.edit2,'String'))
if f1==1
[yd,yn,opt2] = denoise(x,h);
elseif f1==2
 [yd,yn,opt2] = denoise(x,h,1);
else
[yd,yn,opt2] = denoise(x,h);
end



axes(handles.axes4)
plot(x);hold on;plot(yd,'r');hold off;
legend('Noisy Signal','Filter Signal')
set(handles.axes4,'XMinorTick','on')
grid on

uiwait(msgbox('Denoising Done','Denoising','modal'));

else
File=get(handles.edit1,'String');   
x1=wavread(File);
x2=x1(1:65536);
f2 = str2double(get(handles.edit3,'String'))
if f2==4
h = daubcqf(4); 
elseif f2==6
h = daubcqf(6); 
elseif f2==8;
h = daubcqf(8); 
elseif f2==10;
h = daubcqf(10); 
else
h = daubcqf(4); 
end
f1 = str2double(get(handles.edit2,'String'))
if f1==1
[yd,yn,opt2] = denoise(x2,h);
elseif f1==2
 [yd,yn,opt2] = denoise(x2,h,1);
else
[yd,yn,opt2] = denoise(x2,h);
end



axes(handles.axes4)
plot(yd)
set(handles.axes4,'XMinorTick','on')
grid on

uiwait(msgbox('Denoising Done','Denoising','modal'));
end

% --------------------------------------------------------------------
function varargout = File_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Close_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Calculate_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_6_Callback(h, eventdata, handles, varargin)
File=get(handles.edit1,'String');   
x1=wavread(File);
x2=x1(1:65536);
h = daubcqf(6); 
[yd,yn,opt2] = denoise(x2,h,1);
axes(handles.axes2)
plot(yd);
set(handles.axes2,'XMinorTick','on')
grid on

uiwait(msgbox('Denoising Done','Denoising','modal'));




% --------------------------------------------------------------------
function varargout = Help_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_8_Callback(h, eventdata, handles, varargin)
[data map]=imread('sheraz.bmp');
uiwait(msgbox('This Software is developed by Sheraz Khan','About Developer','custom',data,map));


% --------------------------------------------------------------------
function varargout = LoadFile_Callback(h, eventdata, handles, varargin)

function varargout = Open_Callback(h, eventdata, handles, varargin)
% Use UIGETFILE to allow for the selection of a custom address book.
[filename, pathname] = uigetfile( ...
    {'*.wav', 'All wav-Files (*.wav)'; ...
        '*.*','All Files (*.*)'}, ...
    'Select Address Book');
% If "Cancel" is selected then return

    File = fullfile(pathname,filename);
x1=wavread('File');
x2=x1(1:65536);
   h = daubcqf(6); 
   [yd,yn,opt2] = denoise(x2,h,1,[]);
axes(handles.axes1)
plot(x2);hold on;plot(yd,'r');
set(handles.axes1,'XMinorTick','on')
grid on




% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = popupmenu1_Callback(h, eventdata, handles, varargin)

val = get(h,'Value');
switch val
case 1
set(handles.edit2,'String','1')

case 2
set(handles.edit2,'String','2')

end
% --------------------------------------------------------------------
function varargout = popupmenu2_Callback(h, eventdata, handles, varargin)

val = get(h,'Value');
switch val
case 1
set(handles.edit3,'String','4')

case 2
set(handles.edit3,'String','6')
case 3
set(handles.edit3,'String','8')
case 4
set(handles.edit3,'String','10')

end






% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit3_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------



% --------------------------------------------------------------------
function varargout = popupmenu3_Callback(h, eventdata, handles, varargin)
val = get(h,'Value');
switch val
case 1
set(handles.pushbutton4,'Enable','on')
set(handles.text10,'Visible','off')
set(handles.popupmenu4,'Visible','off')
set(handles.text11,'Visible','off')
set(handles.popupmenu5,'Visible','off')
set(handles.radiobutton3,'Visible','on')
set(handles.radiobutton4,'Visible','on')
set(handles.frame2,'Visible','on')
set(handles.text12,'Visible','on')
case 2
set(handles.pushbutton4,'Enable','off')
set(handles.radiobutton3,'Visible','off')
set(handles.radiobutton4,'Visible','off')
set(handles.frame2,'Visible','off')
set(handles.text12,'Visible','off')
set(handles.text10,'Visible','on')
set(handles.popupmenu4,'Visible','on')
set(handles.text11,'Visible','on')
set(handles.popupmenu5,'Visible','on')
end

% --------------------------------------------------------------------
function varargout = popupmenu4_Callback(h, eventdata, handles, varargin)

val = get(h,'Value');
switch val
case 1
set(handles.edit1,'String','Bumps')
signal= 'Bumps';
case 2
set(handles.edit1,'String','Blocks')
signal= 'Blocks';
case 3
set(handles.edit1,'String','Doppler')
signal= 'Doppler';
case 4
set(handles.edit1,'String','Ramp')
signal= 'Ramp';
case 5
set(handles.edit1,'String','Cusp')
signal= 'Cusp';
case 6
set(handles.edit1,'String','Sing')
signal= 'Sing';
case 7
set(handles.edit1,'String','HiSine')
signal= 'HiSine';
case 8
set(handles.edit1,'String','LoSine')
signal= 'LoSine';
case 9
set(handles.edit1,'String','LinChirp')
signal= 'LinChirp';
case 10
set(handles.edit1,'String','TwoChirp')
signal= 'TwoChirp';
case 11
set(handles.edit1,'String','QuadChirp')
signal= 'QuadChirp';
case 12
set(handles.edit1,'String','MishMash')
signal= 'MishMash';
case 13
set(handles.edit1,'String','Werner Sorrows')
signal= 'Werner Sorrows';
case 14
set(handles.edit1,'String','Leopold')
signal= 'Leopold';
case 15
set(handles.edit1,'String','HeaviSine')
signal= 'HeaviSine';
end
f2 = str2double(get(handles.edit4,'String'));
if f2==2
h = 2; 
elseif f2==4
h = 4; 
elseif f2==6;
h = 6; 
elseif f2==8;
h = 8 ; 
elseif f2==10;
h = 10; 
else
h = 2; 
end  

[s,N] = makesig(signal);
   n = randn(1,N);
x = s + n/h; 
    axes(handles.axes1)
plot(x);hold on;plot(s,'r');hold off;
legend('Noisy Signal','Pure Signal')
set(handles.axes1,'XMinorTick','on')
grid on



% --------------------------------------------------------------------
function varargout = popupmenu5_Callback(h, eventdata, handles, varargin)
val = get(h,'Value');
switch val
case 1
set(handles.edit4,'String','2')
case 2
set(handles.edit4,'String','4')
case 3
set(handles.edit4,'String','6')
case 4
set(handles.edit4,'String','8')
case 5
set(handles.edit4,'String','10')
end



% --------------------------------------------------------------------
function varargout = edit4_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit5_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = radiobutton3_Callback(h, eventdata, handles, varargin)
set(handles.radiobutton4,'Value',0)
set(handles.axes1,'Visible','on')
set(handles.axes7,'Visible','off')
x2=handles.X1;
fs=handles.X2
axes(handles.axes1)
plot(x2);
set(handles.axes1,'XMinorTick','on')
grid on



% --------------------------------------------------------------------
function varargout = radiobutton4_Callback(h, eventdata, handles, varargin)

set(handles.radiobutton3,'Value',0)
set(handles.axes1,'Visible','off')
set(handles.axes7,'Visible','on')
x2=handles.X1;
fs=handles.X2
fft1=fft(x2,65536);
w=(0:32767)/32768*(fs/2);
axes(handles.axes7)
plot(w,abs(fft1(1:32768)));
set(handles.axes7,'XMinorTick','on')
grid on


% --------------------------------------------------------------------
function varargout = Signal_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Length_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_12_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Type_Callback(h, eventdata, handles, varargin)

