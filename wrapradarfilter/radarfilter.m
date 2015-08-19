function [ret_observations, ao] = radarfilter(observations)

ret_observations = observations;
cfg_filename 					= observations('cfg_filename');
additional_output_filename 		= observations('additional_output_filename');
radar_azel_r1_m 				= observations('radar_azel_r1_m');
radar_azel_r2_m 				= observations('radar_azel_r2_m');
radar_azel_alpha_rad 			= observations('radar_azel_alpha_rad');
radar_azel_gamma_rad 			= observations('radar_azel_gamma_rad');
radar_beam_FWHM0_rad 			= observations('radar_beam_FWHM0_rad');
radar_beam_FWHM1_rad 			= observations('radar_beam_FWHM1_rad');
radar_t 						= observations('radar_t');
radar_dt 						= observations('radar_dt');
radar_SNR 						= observations('radar_SNR');

[mypathstr,myname,myext] = fileparts(mfilename('fullpath'));
oldFolder = cd(mypathstr);

addpath(strcat(mypathstr, '..\..\additional_output'));

fileID = fopen('./input.z', 'w');
fprintf(fileID, '%s %s\n', 'cfg_filename', cfg_filename);
fprintf(fileID, '%s %s\n', 'additional_output_filename', additional_output_filename);

fprintf(fileID, '%-30s %-15i', 'radar_azel_r1_m', length(radar_azel_r1_m));
for j=1:length(radar_azel_r1_m)
    fprintf(fileID, ' %-15.3e', radar_azel_r1_m(j));
end
fprintf(fileID, '\n');

fprintf(fileID, '%-30s %-15i', 'radar_azel_r2_m', length(radar_azel_r2_m));
for j=1:length(radar_azel_r2_m)
    fprintf(fileID, ' %-15.3e', radar_azel_r2_m(j));
end
fprintf(fileID, '\n');

fprintf(fileID, '%-30s %-15i', 'radar_azel_alpha_rad', length(radar_azel_alpha_rad));
for j=1:length(radar_azel_alpha_rad)
    fprintf(fileID, ' %-15.3e', radar_azel_alpha_rad(j));
end
fprintf(fileID, '\n');

fprintf(fileID, '%-30s %-15i', 'radar_azel_gamma_rad', length(radar_azel_gamma_rad));
for j=1:length(radar_azel_gamma_rad)
    fprintf(fileID, ' %-15.3e', radar_azel_gamma_rad(j));
end
fprintf(fileID, '\n');

fprintf(fileID, '%-30s %-15i', 'radar_beam_FWHM0_rad', length(radar_beam_FWHM0_rad));
for j=1:length(radar_beam_FWHM0_rad)
    fprintf(fileID, ' %-15.3e', radar_beam_FWHM0_rad(j));
end
fprintf(fileID, '\n');

fprintf(fileID, '%-30s %-15i', 'radar_beam_FWHM1_rad', length(radar_beam_FWHM1_rad));
for j=1:length(radar_beam_FWHM1_rad)
    fprintf(fileID, ' %-15.3e', radar_beam_FWHM1_rad(j));
end
fprintf(fileID, '\n');

fprintf(fileID, '%-30s %-15i', 'radar_t', length(radar_t));
for j=1:length(radar_t)
    fprintf(fileID, ' %-15.3e', radar_t(j));
end
fprintf(fileID, '\n');

fprintf(fileID, '%-30s %-15i', 'radar_dt', length(radar_dt));
for j=1:length(radar_dt)
    fprintf(fileID, ' %-15.3e', radar_dt(j));
end
fprintf(fileID, '\n');

fprintf(fileID, '%-30s %-15i', 'radar_SNR', length(radar_SNR));
for j=1:length(radar_SNR)
    fprintf(fileID, ' %-15.3e', radar_SNR(j));
end
fprintf(fileID, '\n');

fclose(fileID);

system('wrapradarfilter_matlab.bat');

fileID = fopen('./output.z');

fscanf(fileID,'%*s %i', 1);
radar_dBZ = fscanf(fileID,'%f', [length(radar_azel_r1_m)]);
fscanf(fileID,'%*s %i', 1);
radar_radial_vel_ms = fscanf(fileID,'%f', [length(radar_azel_r1_m)]);
fscanf(fileID,'%*s %i', 1);
radar_radial_vel_width_ms = fscanf(fileID,'%f', [length(radar_azel_r1_m)]);
fscanf(fileID,'%*s %i', 1);
radial_vel_shape = fscanf(fileID,'%f', [8*length(radar_azel_r1_m)]);

fclose(fileID);

ret_observations('radar_dBZ') 					= radar_dBZ;
ret_observations('radar_radial_vel_ms') 		= radar_radial_vel_ms;
ret_observations('radar_radial_vel_width_ms') 	= radar_radial_vel_width_ms;
ret_observations('radial_vel_shape') 			= radial_vel_shape;

ret_observations('nspectrum') = 8;

ret_observations('radial_vel_shape') = reshape(ret_observations('radial_vel_shape'), [ret_observations('nspectrum'), ret_observations('n')]);
ret_observations('radial_vel_shape') = transpose(ret_observations('radial_vel_shape'));

ao = additional_output(observations('additional_output_filename'));

cd(oldFolder);
end