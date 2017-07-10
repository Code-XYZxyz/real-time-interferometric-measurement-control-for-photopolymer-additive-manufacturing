function m = icm_init_mem_file(cp, clearFile)

filename = fullfile(pwd, 'icm_comm.dat');

numFitter     = cp.numFitter;
W             = cp.imW;
H             = cp.imH;
numImageCache = cp.maxNumCache;

frame = zeros(H, W);
frameIdx = 0;
frameTime = 0.0;
tic = 0;

% global status, only use first one of the array
%  (1): measurement status. 1: running, 0 stopped.
%  (2): lastWriteFrameIdx
%  (3): workerStatus. 0: not running, 1: started, 2: targetReached
status = [1 0 0]';

if clearFile && exist(filename, 'file')
    disp(['remove file ', filename])
    delete(filename)
end

% Create the communications file if it is not already there.
if ~exist(filename, 'file')
    disp(['init file ', filename])
    [f, msg] = fopen(filename, 'wb');
    if f ~= -1
        for i = 1:numImageCache
            fwrite(f, frame, 'uint8');
            fwrite(f, frameIdx, 'uint32');
            fwrite(f, frameTime, 'double');
            fwrite(f, 0, 'uint64'); % snapTic
            fwrite(f, 0, 'uint8'); % uvIris
            fwrite(f, 0, 'uint8'); % uvStatus
            fwrite(f, 0, 'double'); % avgTotalHeight
            fwrite(f, 0, 'double'); % avgTotalPhase
            fwrite(f, status, 'int32');
        end
        fclose(f);
    else
        error('MATLAB:demo:send:cannotOpenFile', ...
              'Cannot open file "%s": %s.', filename, msg);
    end
end
 
% Memory map the file.
disp(['creating memmapfile for', filename])
m = memmapfile(filename,...
    'Writable', true,...
    'Format', {...
        'uint8', size(frame), 'frame';...
        'uint32', size(frameIdx), 'frameIdx';...
        'double', size(frameTime), 'frameTime';...
        'uint64', [1 1], 'snapTic';...
        'uint8', [1 1], 'uvIris';...
        'uint8', [1 1], 'uvStatus';...
        'double', [1, 1], 'avgTotalHeight';...
        'double', [1, 1], 'avgTotalPhase';...
        'int32', size(status), 'status';...
    },...
    'Repeat', numImageCache...
    );

fprintf('init status = %d\n', m.Data(1).status);

return;
