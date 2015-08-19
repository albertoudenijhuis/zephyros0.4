addpath('..\..\..\wrapwindfield');


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
	
    cfg_filename =  '..\input_files\simulation\instruments\simulation_tara.cfg;';
    cfg_filename = strcat(cfg_filename,'..\input_files\simulation\scatterers\homogeneous.cfg;');
    cfg_filename = strcat(cfg_filename,'..\input_files\simulation\wind\',myplot,'.cfg;');
    
    additional_output_filename = '.\additional_output.zout';

    
    for k=1:3,
        if k == 1,
           plane = 'xy';
           xvec = [-10.e3:1.e3:10.e3];
           yvec = [-10.e3:1.e3:10.e3];
           zvec = [0];
           tvec = [0];
        end
        if k == 2,
           plane = 'xz';
           xvec = [-10.e3:1.e3:10.e3];
           yvec = [0];
           zvec = [0.:1.e3:10.e3];
           tvec = [0];
        end
        if k == 3,
           plane = 'yz';
           xvec = [-10.e3:1.e3:10.e3];
           yvec = [0];
           zvec = [0.:1.e3:10.e3];
           tvec = [0];
        end

        [x,y,z,t] = ndgrid(xvec,yvec,zvec,tvec);
        x = reshape(x,[],1);
        y = reshape(y,[],1);
        z = reshape(z,[],1);
        t = reshape(t,[],1);
        
        [u,v,w] = windfield(cfg_filename, additional_output_filename, x, y, z, t);
        
       f = figure('visible','off');
       if plane == 'xy',
           xlabel('x');
           ylabel('y');
           lh = quiver(x,y,u,v);
       end
       if plane == 'xz',
           xlabel('x');
           ylabel('z');
           lh = quiver(x,z,u,w);
       end
       if plane == 'yz',
           xlabel('y');
           ylabel('z');
           lh = quiver(x,z,u,w);
       end
       
       myname = strcat('plots\windfield_',myplot,'_',plane,'.png')
       print(f, '-r80', '-dpng', myname);   

    end    
end
    
 
