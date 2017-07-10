function handles = updateParameters(handles)
handles.cp.FPS = str2double(get(handles.ed_FPS,'String'));
handles.cp.MeasPeriodSamples = str2double(get(handles.ed_MeasPeriodSamples,'String'));
handles.cp.SamplesNumB4Measure = str2double(get(handles.ed_SamplesNumBeforeMeasure,'String'));
handles.cp.MovingHorizonL = str2double(get(handles.ed_MovingHorizonL,'String'));
handles.cp.HalfLife = str2double(get(handles.ed_HalfLife,'String'));
handles.cp.GOF_rSquare = str2double(get(handles.ed_GOF_rSquare,'String'));
handles.cp.Wavelength = str2double(get(handles.ed_Wavelength,'String'));
handles.cp.n_L = str2double(get(handles.ed_n_L,'String'));
handles.cp.n_m = str2double(get(handles.ed_n_m,'String'));
handles.cp.f_max = str2double(get(handles.ed_f_max,'String'));
handles.cp.f_diff_max = str2double(get(handles.ed_f_diff_max,'String'));

% In batch job, pool of 1 main worker + N fitter
handles.cp.numFitter = str2double(get(handles.ed_NumFitters,'String'));

handles.cp.uvIris = str2double(get(handles.ed_uvIris,'String'));
% mode 0: None, 1: height, 2: time
handles.cp.targetMode = 0;
if get(handles.rb_TargetCuredHeight, 'Value') == 1.0
    handles.cp.targetMode = 1;
    handles.cp.targetCuredHeight = str2double(get(handles.ed_TargetCuredHeight, 'String'));
elseif get(handles.rb_TargetExpTime, 'Value') == 1.0
    handles.cp.targetMode = 2;
    handles.cp.targetExpTime = str2double(get(handles.ed_TargetExpTime, 'String'));
end

end
