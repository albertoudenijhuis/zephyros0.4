#!/usr/bin/env python
import re
import numpy as np

def ao_dct(file1):
    dct = {}
    for linestr in open(file1,'r').readlines():
        if linestr[:2] == '!!':
            rawdata = linestr[2:]
            rawdata = rawdata.split()
            
            #varname = linestr[2:32].strip()
            #rawdata = re.findall('(\-?[0-9]{1,10}\.?[0-9]*E?e?\+?\-?[0-9]*|nan|NaN|inf)' , linestr[2:] )
            varname = rawdata[0].strip()
            ndim    = int(rawdata[1])
            dims    = list(map(int,rawdata[2:2+ndim]))
            rawdata2 = list(map(float, rawdata[2+ndim:]))
            #check shape
            if len(rawdata2) == np.product(dims):
                dct[varname] = np.array(rawdata2).reshape(dims, order='C')
            else:
                print("Problem with reading out additional output. Variable ", varname, " has length ", len(rawdata2), " and dimensions ", dims)
    return dct


#write additional output or measurement file
#in development
def write_additional_output(dct, filename):
    #txtzeros = build_string(np.zeros(data['dBZ_hh_err'].shape))
    #txtones = build_string(np.ones(data['dBZ_hh_err'].shape))
    
    preferred_order = [
        'center_x',
        'center_y',
        'center_z',
        'center_t',
        'enu_radar_location_x',
        'enu_radar_location_y',
        'enu_radar_location_z',
        'enu_radar_time_t',
        'azel_r1_m',
        'azel_r2_m',
        'azel_alpha_rad',
        'azel_gamma_rad',
        'beam_FWHM0_rad',
        'beam_FWHM1_rad',
        'dBZ_hh',
        'dBZ_hh_err',
        'Doppler_velocity_hh_ms',
        'Doppler_velocity_hh_ms_err',
        'Doppler_spectral_width_hh_ms',
        'Doppler_spectral_width_hh_ms_err',
    ]
    
    
    file_txt = ""

    #item from preferred order
    for item in preferred_order:
        if item in dct.keys():
            file_txt += "!! "+item+"           "+build_string(dct[item])+"\n"

    #~ #other items in sorted order
    #~ for item in sorted(dct.keys()):
        #~ if not (item in preferred_order):
            #~ file_txt += "!! "+item+"           "+build_string(dct[item])+"\n"   
    
    f = open(filename, "w")
    f.write(file_txt)
    f.close()





def build_string(var):
    txt = "{:<15}".format(len(var.shape))   
    for i in range(len(var.shape)):
        txt += "{:<15}".format(var.shape[i])     
    var2 = np.reshape(var, [np.product(var.shape)], order='C')
    for i in range(len(var2)):
        txt   += "{:<15.3e}".format(var2[i])

    return txt




#solution to find back identifiers
def search_identifier(file1, section, subsection, identifier, after_linenr = 0):
    linenr = 0
    line_section    = ""
    line_subsection = ""
    line_identifier = ""
    
    for linestr in open(file1,'r').readlines():
        linenr += 1

        #ignore comments
        if ((len(linestr) > 0) and (linestr[0] == "#")):
            continue

        if len(linestr) == 0:
            continue
            
        line_splitted = linestr.split()
        
        if (len(line_splitted) > 1):            
            if (line_splitted[0] == 'section'):
                line_section = line_splitted[1]
            elif (line_splitted[0] == 'subsection'):
                line_subsection = line_splitted[1]
                
        #ignore some lines
        if (linenr <= after_linenr):
            continue
                
        if (len(line_splitted) > 1):                           
            line_identifier = line_splitted[0]
            if (
                (section == line_section) &
                (subsection == line_subsection) &
                (identifier == line_identifier)
                ):
                return linenr, linestr

    return False, False
