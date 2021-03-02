function [ret_observations, ao] = radarfilter(observations)

ret_observations = observations;
cfg_filename 					= observations('cfg_filename');
additional_output_filename 		= observations('additional_output_filename');
o_enu_x0 						= observations('o_enu_x0');
o_enu_y0 						= observations('o_enu_y0');
o_enu_z0 						= observations('o_enu_z0');
o_azel_r1_m 					= observations('o_azel_r1_m');
o_azel_r2_m 					= observations('o_azel_r2_m');
o_azel_alpha_rad 				= observations('o_azel_alpha_rad');
o_azel_gamma_rad 				= observations('o_azel_gamma_rad');
o_beam_FWHM0_rad 				= observations('o_beam_FWHM0_rad');
o_beam_FWHM1_rad 				= observations('o_beam_FWHM1_rad');
o_t				 				= observations('o_t');
o_dt			 				= observations('o_dt');
o_refl			 				= observations('o_refl');
o_srefl			 				= observations('o_srefl');
o_vr			 				= observations('o_vr');
o_svr			 				= observations('o_svr');
o_spectralwidth			 		= observations('o_spectralwidth');
o_sspectralwidth			 	= observations('o_sspectralwidth');

[mypathstr,myname,myext] = fileparts(mfilename('fullpath'));
oldFolder = cd(mypathstr);

addpath(strcat(mypathstr, '..\..\additional_output'));

fileID = fopen('./input.z', 'w');
fprintf(fileID, '%s %s\r\n', 'cfg_filename', cfg_filename);
fprintf(fileID, '%s %s\r\n', 'additional_output_filename', additional_output_filename);

fprintf(fileID, '%-30s %-15i', 'o_enu_x0', length(o_enu_x0));
for j=1:length(o_enu_x0)
    fprintf(fileID, ' %-15.3e', o_enu_x0(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_enu_y0', length(o_enu_y0));
for j=1:length(o_enu_y0)
    fprintf(fileID, ' %-15.3e', o_enu_y0(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_enu_z0', length(o_enu_z0));
for j=1:length(o_enu_z0)
    fprintf(fileID, ' %-15.3e', o_enu_z0(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_azel_r1_m', length(o_azel_r1_m));
for j=1:length(o_azel_r1_m)
    fprintf(fileID, ' %-15.3e', o_azel_r1_m(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_azel_r2_m', length(o_azel_r2_m));
for j=1:length(o_azel_r2_m)
    fprintf(fileID, ' %-15.3e', o_azel_r2_m(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_azel_alpha_rad', length(o_azel_alpha_rad));
for j=1:length(o_azel_alpha_rad)
    fprintf(fileID, ' %-15.3e', o_azel_alpha_rad(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_azel_gamma_rad', length(o_azel_gamma_rad));
for j=1:length(o_azel_gamma_rad)
    fprintf(fileID, ' %-15.3e', o_azel_gamma_rad(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_beam_FWHM0_rad', length(o_beam_FWHM0_rad));
for j=1:length(o_beam_FWHM0_rad)
    fprintf(fileID, ' %-15.3e', o_beam_FWHM0_rad(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_beam_FWHM1_rad', length(o_beam_FWHM1_rad));
for j=1:length(o_beam_FWHM1_rad)
    fprintf(fileID, ' %-15.3e', o_beam_FWHM1_rad(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_t', length(o_t));
for j=1:length(o_t)
    fprintf(fileID, ' %-15.3e', o_t(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_dt', length(o_dt));
for j=1:length(o_dt)
    fprintf(fileID, ' %-15.3e', o_dt(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_refl', length(o_refl));
for j=1:length(o_refl)
    fprintf(fileID, ' %-15.3e', o_refl(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_srefl', length(o_srefl));
for j=1:length(o_srefl)
    fprintf(fileID, ' %-15.3e', o_srefl(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_vr', length(o_vr));
for j=1:length(o_vr)
    fprintf(fileID, ' %-15.3e', o_vr(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_svr', length(o_svr));
for j=1:length(o_svr)
    fprintf(fileID, ' %-15.3e', o_svr(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_spectralwidth', length(o_spectralwidth));
for j=1:length(o_spectralwidth)
    fprintf(fileID, ' %-15.3e', o_spectralwidth(j));
end
fprintf(fileID, '\r\n');

fprintf(fileID, '%-30s %-15i', 'o_sspectralwidth', length(o_sspectralwidth));
for j=1:length(o_sspectralwidth)
    fprintf(fileID, ' %-15.3e', o_sspectralwidth(j));
end
fprintf(fileID, '\r\n');

fclose(fileID);

system('wrapwindvectors_matlab.bat');

fileID = fopen('./output.z');

thisline = fgetl(fileID); thisdata = strsplit(strtrim(thisline));
windvector_u = reshape(cellfun(@(x) str2num(x), thisdata(3:length(thisdata))), [], 1);
thisline = fgetl(fileID); thisdata = strsplit(strtrim(thisline));
windvector_su = reshape(cellfun(@(x) str2num(x), thisdata(3:length(thisdata))), [], 1);
thisline = fgetl(fileID); thisdata = strsplit(strtrim(thisline));
windvector_v = reshape(cellfun(@(x) str2num(x), thisdata(3:length(thisdata))), [], 1);
thisline = fgetl(fileID); thisdata = strsplit(strtrim(thisline));
windvector_sv = reshape(cellfun(@(x) str2num(x), thisdata(3:length(thisdata))), [], 1);
thisline = fgetl(fileID); thisdata = strsplit(strtrim(thisline));
windvector_w = reshape(cellfun(@(x) str2num(x), thisdata(3:length(thisdata))), [], 1);
thisline = fgetl(fileID); thisdata = strsplit(strtrim(thisline));
windvector_sw = reshape(cellfun(@(x) str2num(x), thisdata(3:length(thisdata))), [], 1);

fclose(fileID);

ret_observations('windvector_u') 				= windvector_u;
ret_observations('windvector_su') 				= windvector_su;
ret_observations('windvector_v') 				= windvector_v;
ret_observations('windvector_sv') 				= windvector_sv;
ret_observations('windvector_w') 				= windvector_w;
ret_observations('windvector_sw') 				= windvector_sw;

ao = additional_output(observations('additional_output_filename'))

cd(oldFolder);
end
