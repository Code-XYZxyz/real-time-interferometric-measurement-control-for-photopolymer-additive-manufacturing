% --- Executes on button press in pb_DefaultMeasParameters.
function pb_DefaultMeasParameters_Callback(hObject, eventdata, handles)
% handles = guidata(hObject);
set(handles.ed_FPS, 'string', '30');
set(handles.ed_MeasPeriodSamples, 'string', '10');

% used 20, then found 5 could detect curing more timely
set(handles.ed_SamplesNumBeforeMeasure,'string','0');

set(handles.ed_MovingHorizonL, 'string', '32');
% set(handles.ed_MovingHorizonL, 'string', '48');

% half life use 10 or 
% ..."handles.ed_MeasPeriodSamples" if which is more than 10
% set(handles.ed_HalfLife, 'string', num2str(max(10,str2num(get(handles.ed_MeasPeriodSamples,'String')))));
set(handles.ed_HalfLife, 'string','10');

set(handles.ed_GOF_rSquare,'string','0.95'); % 0.85 could detect threshold better
set(handles.ed_Wavelength, 'string', '0.532');
set(handles.ed_n_L, 'string', '1.4723');

% initial guess: 1.5002
...03/15/2016: withoud dark curing: 1.4946, this could actually be larger once dark curing included
% set(handles.ed_n_m, 'string', '1.4946');
set(handles.ed_n_m, 'string', '1.4945');

% Stopwatch control: time or height
set(handles.rb_TargetCuredHeight,'Value', 1); % Target height
set(handles.ed_TargetCuredHeight,'string', 40); % Target height
set(handles.rb_TargetExpTime,'Value', 0); % Target Time
set(handles.ed_TargetExpTime,'string', 40); % Target Time


% to detect outlier frequency
% maximum frequency (Hz):1.5
% set(handles.ed_f_max,'string','1.5');
set(handles.ed_f_max,'string','15');
% maximum frequency difference
set(handles.ed_f_diff_max,'string','0.6');

% In batch job, pool of 1 main worker + N fitter
handles.cp.numFitter = set(handles.ed_NumFitters,'String','2');

% uvIris
set(handles.ed_uvIris,'String','5');

guidata(hObject, handles);
