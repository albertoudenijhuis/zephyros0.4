addpath('..\..\..\wrapradarfilter');

for j=1:5,
    if j == 1,
        myplot = 'vector';
    end
    if j == 2,
        myplot = 'grid';
    end
    if j == 3,
        myplot = 'wave';
    end
    if j == 4,
        myplot = 'rankine_vortex';
    end
    if j == 5,
        myplot = 'lamb_oseen_vortex';
    end

    observations = containers.Map;

    observations('cfg_filename')	=  '..\input_files\simulation\instruments\simulation_tara.cfg;';
    observations('cfg_filename')	 = strcat(observations('cfg_filename')	,'..\input_files\simulation\scatterers\homogeneous.cfg;');
    observations('cfg_filename')	 = strcat(observations('cfg_filename')	,'..\input_files\simulation\wind\',myplot,'.cfg;');

    observations('additional_output_filename') = '.\additional_output.zout';

    dr = 30.;
    elevation_angle = 45.;
    observations('radar_azel_r1_m') 		= 0:dr:15.e3;
    observations('n')  						= length(observations('radar_azel_r1_m'));
    observations('radar_azel_r2_m') 		= observations('radar_azel_r1_m') + dr;
    observations('radar_azel_alpha_rad')    = zeros(observations('n'),1) + (pi/180.) * 45.;
    observations('radar_azel_gamma_rad')    = zeros(observations('n'),1) + (pi/180.) * elevation_angle;
    observations('radar_beam_FWHM0_rad')  	= zeros(observations('n'),1) + (pi/180.) * 10.;
    observations('radar_beam_FWHM1_rad')  	= zeros(observations('n'),1) + (pi/180.) * 10.;
    observations('radar_t')          		= zeros(observations('n'),1);
    observations('radar_dt')          		= ones(observations('n'),1) ;
    observations('radar_SNR')          		= zeros(observations('n'),1) + 1000.;

    [observations, ao] = radarfilter(observations);


    %obtain velocity-values of spectrum
    mincdfP = 1.e-50;
    maxcdfP = 1. - 1.e-50;    
    cdfvaluesa = linspace(mincdfP, maxcdfP, observations('nspectrum') + 1);
    cdfvalues = (0.5 * cdfvaluesa(1:end-1)) + (0.5 * cdfvaluesa(2:end));
    fstat_values_xx = norminv(cdfvalues, 0., 1.);
    fstat_values_dx = norminv(cdfvaluesa(2:end), 0.,1.) - norminv(cdfvaluesa(1:end-1), 0., 1.);

    velocities = zeros(observations('n') , observations('nspectrum'));
    ranges = zeros(observations('n') , observations('nspectrum'));

    o_range = observations('radar_azel_r1_m') + dr/2.;
    rad_vel                 = observations('radar_radial_vel_ms');
    radial_vel_width_ms     = observations('radar_radial_vel_width_ms');
    rad_shape               = observations('radial_vel_shape');

    for i_m=1:observations('n'),
        for i_s=1:observations('nspectrum'),
            %i = (i_m - 1) * nspectrum + i_s;
          ranges(i_m, i_s)         = o_range(i_m);
          velocities(i_m, i_s)     = rad_vel(i_m) + radial_vel_width_ms(i_m) * fstat_values_xx(i_s);
          rad_shape(i_m, i_s)      = rad_shape(i_m, i_s) /  fstat_values_dx(i_s);
    end
    end

    f = figure('visible','off');
    [C,h] = contourf(velocities,ranges,rad_shape);
    set(h,'LineColor','none')

    xlabel('velocity [m/s]');
    ylabel('range [m]');
    colorbar;
    myname = strcat('radarfilter_staring_',myplot);
    print(f, '-r80', '-dpng', strcat('plots\', myname, '.png'));  
	save(strcat('data\', myname, '.mat'), 'observations', 'ao');
end
