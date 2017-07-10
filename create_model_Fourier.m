function [ft, opts] = create_model_Fourier(dtrain, wgt)
% y = a0 + a1*cos(w*x) + b1*sin(w*x)
% Coefficients = [a0;a1;b1;w];
ft = fittype('fourier1');

opts = fitoptions(ft);
opts.Weights = wgt';

opts.Display = 'Off';

%% note: setting reasonable limit is very important!
% % set limits of the coefficients "a0,a1,b1,w"
wMax = 1.0/15;
opts.Lower = [0 -255 -255 0];
opts.Upper = [255 255 255 wMax];
% opts.StartPoint = [0 0 0 0];

% TODO:
% 1. Use rsquare to check fit goodness
% 2. Use previous good parameter as starting point of next fit
% 3. Adaptive half life and moving horizon if rsquare too low

end