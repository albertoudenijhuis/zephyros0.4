addpath('..\..\..\wrapwindvectors');

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

    myname = strcat('radarfilter_scanning_',myplot);
	load(strcat('..\..\wrapradarfilter\matlab\data\', myname, '.mat'), 'observations', 'ao');
	
    observations('cfg_filename')	=  '..\input_files\retrieval\windvectors\windvectors_lwm.cfg';
    observations('additional_output_filename') = '.\additional_output.zout';

	observations('o_enu_x0')				= zeros(observations('n'),1);
	observations('o_enu_y0')				= zeros(observations('n'),1);
	observations('o_enu_z0')				= zeros(observations('n'),1);
	observations('o_azel_r1_m')				= observations('radar_azel_r1_m');
	observations('o_azel_r2_m')				= observations('radar_azel_r2_m');
	observations('o_azel_alpha_rad')		= observations('radar_azel_alpha_rad');
	observations('o_azel_gamma_rad')		= observations('radar_azel_gamma_rad');
	observations('o_beam_FWHM0_rad')		= observations('radar_beam_FWHM0_rad');
	observations('o_beam_FWHM1_rad')		= observations('radar_beam_FWHM1_rad');
	observations('o_t')						= observations('radar_t');
	observations('o_dt')					= observations('radar_dt');
	observations('o_refl')					= zeros(observations('n'),1);
	observations('o_srefl')					= zeros(observations('n'),1) + 1.;
	observations('o_vr')					= observations('radar_radial_vel_ms');
	observations('o_svr')					= zeros(observations('n'),1) + 1.;
	observations('o_spectralwidth')			= zeros(observations('n'),1) + 1.;
	observations('o_sspectralwidth')		= zeros(observations('n'),1) + 1.;
	
    [observations, ao] = windvectors(observations);

    myname = strcat('windvectors_scanning_',myplot);

    f = figure('visible','off');
    hold on    
    xlin=linspace(min(ao('x')),max(ao('x')),20);
    ylin=linspace(min(ao('y')),max(ao('y')),20);
    [X,Y]=meshgrid(xlin,ylin);
    U=griddata(ao('x'),ao('y'),observations('windvector_u'),X,Y,'cubic');
    V=griddata(ao('x'),ao('y'),observations('windvector_v'),X,Y,'cubic');
    Z=griddata(ao('x'),ao('y'),ao('lwm_vr'),X,Y,'cubic');

    set(gca,'color','none')
    [C,h] = contourf(X,Y,Z, 50);
    set(h,'LineColor','none')
    
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
    
    maxval = max(abs( ao('lwm_vr')));
    caxis([-maxval, maxval]);
    lh = quiver(X,Y,U,V);
    set(lh,'linewidth',2);
    set(lh,'color',[1,1,1]);
     
    xlabel('x [m]');
    ylabel('y [m]');
    print(f, '-r80', '-dpng', strcat('plots\', myname, '.png'));  
    hold off % reset hold state

end
