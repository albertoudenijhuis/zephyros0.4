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
    di = 10;
    elevation_angle = 0.;
    observations('vec_radar_azel_r1_m') 		= 1:di*dr:15.e3;
    observations('vec_radar_azel_alpha_rad')    = (pi/180.) * (0:5:355);

    [observations('radar_azel_r1_m'),observations('radar_azel_alpha_rad')] = ndgrid(observations('vec_radar_azel_r1_m'),observations('vec_radar_azel_alpha_rad'));
    observations('radar_azel_r1_m') = reshape(observations('radar_azel_r1_m'),[],1);
    observations('radar_azel_alpha_rad') = reshape(observations('radar_azel_alpha_rad'),[],1);

    observations('n')  						= length(observations('radar_azel_r1_m'));
    observations('radar_azel_r2_m') 		= observations('radar_azel_r1_m') + dr;
    observations('radar_azel_gamma_rad')    = zeros(observations('n'),1) + (pi/180.) * elevation_angle;
    observations('radar_beam_FWHM0_rad')  	= zeros(observations('n'),1) + (pi/180.) * 10.;
    observations('radar_beam_FWHM1_rad')  	= zeros(observations('n'),1) + (pi/180.) * 10.;
    observations('radar_t')          		= zeros(observations('n'),1);
    observations('radar_dt')          		= ones(observations('n'),1) ;
    observations('radar_SNR')          		= zeros(observations('n'),1) + 1000.;

    [observations, ao] = radarfilter(observations);


    f = figure('visible','off');
    xlin=linspace(min(ao('x')),max(ao('x')),100);
    ylin=linspace(min(ao('y')),max(ao('y')),100);
    [X,Y]=meshgrid(xlin,ylin);
    Z=griddata(ao('x'),ao('y'),observations('radar_radial_vel_ms'),X,Y,'cubic');
    
    [C,h] = contourf(X,Y,Z, 50);
    set(h,'LineColor','none')

    xlabel('x [m]');
    ylabel('y [m]');
    colorbar;
    b = [0 0 1];       %# start
    w = [.9 .9 .9];    %# middle
    r = [1 0 0];       %# end

    %# colormap of size 64-by-3, ranging from blue -> white -> red
    c1 = zeros(32,3); c2 = zeros(32,3);
    for i=1:3
        c1(:,i) = linspace(b(i), w(i), 32);
        c2(:,i) = linspace(w(i), r(i), 32);
    end
    c = [c1(1:end-1,:);c2];
    colormap(c);
    
    maxval = max(abs(observations('radar_radial_vel_ms')));
    caxis([-maxval, maxval]);
    myname = strcat('radarfilter_scanning_',myplot);
    print(f, '-r80', '-dpng', strcat('plots\', myname, '.png'));  
	save(strcat('data\', myname, '.mat'), 'observations', 'ao');

end
