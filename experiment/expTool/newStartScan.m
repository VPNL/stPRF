function [status, time0] = newStartScan
% This code will trigger the 3T GE scanner at CNI using the E-prime trigger
% cable. --Michael 09/30/2013
try
    s = serial('/dev/tty.usbmodem12341', 'BaudRate', 57600);
    fopen(s);
    fprintf(s, '[t]');
    fclose(s);
catch err
end

if exist('err','var') == 0
    time0  = GetSecs;
    status = 0;
else
    time0 = GetSecs;
    status = 1;
end
end