function [u, v, w] = windfield(cfg_filename, additional_output_filename, x, y, z, t)

[mypathstr,myname,myext] = fileparts(mfilename('fullpath'));
oldFolder = cd(mypathstr);

fileID = fopen('./input.z', 'w');
fprintf(fileID, '%s %s\r\n', 'cfg_filename', cfg_filename);
fprintf(fileID, '%s %s\r\n', 'additional_output_filename', additional_output_filename);

fprintf(fileID, '%-30s %-15i', 'x', length(x));
for j=1:length(x)
    fprintf(fileID, ' %-15.3e', x(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'y', length(y));
for j=1:length(y)
    fprintf(fileID, ' %-15.3e', y(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'z', length(z));
for j=1:length(z)
    fprintf(fileID, ' %-15.3e', z(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 't', length(t));
for j=1:length(t)
    fprintf(fileID, ' %-15.3e', t(j));
end
fprintf(fileID, '\r\n');

fclose(fileID);

system('wrapwindfield_matlab.bat');

fileID = fopen('./output.z');

thisline = fgetl(fileID); thisdata = strsplit(strtrim(thisline));
u = reshape(cellfun(@(x) str2num(x), thisdata(3:length(thisdata))), [], 1);
thisline = fgetl(fileID); thisdata = strsplit(strtrim(thisline));
v = reshape(cellfun(@(x) str2num(x), thisdata(3:length(thisdata))), [], 1);
thisline = fgetl(fileID); thisdata = strsplit(strtrim(thisline));
w = reshape(cellfun(@(x) str2num(x), thisdata(3:length(thisdata))), [], 1);




fclose(fileID);

cd(oldFolder);
end
