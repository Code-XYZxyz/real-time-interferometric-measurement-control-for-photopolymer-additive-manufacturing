function varargout = Real_Time_ICM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%% ----------beginning of xyz inserted code----------
handles.output = hObject; % xyz inserted code
%% ------------end of xyz inserted codes---------------
varargout{1} = handles.output;
