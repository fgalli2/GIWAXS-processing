import os
import datetime
import configparser
import argparse

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors


##############################################################################
class ReadingInputFunctions:
    """Class containing the functions to read the configuation file."""
    def __init__(self):
        """Initialize the class setting the configuration file."""
        parser = argparse.ArgumentParser(
                description="Module to process GIWAXS data")
        parser.add_argument('--config', required=False, type=str,
                            default='configuration_settings',
                            help=('Enter the name of the configuration file '
                                  + '(without .ini extension)'))
        
        config = configparser.ConfigParser()
        args = parser.parse_args()
        config_name = args.config
        config.read(config_name+'.ini')
        self.config = config

    def read_input_files(self):
        section_name = 'Input_files'
        
        data_main_folder = self.config.get(section_name, 'data_main_folder')
        sample_name = self.config.get(section_name, 'sample')
        working_folder = data_main_folder + '\\' + sample_name
        
        calibration_folder = self.config.get(section_name, 'calibration_folder')
        poni_name = self.config.get(section_name, 'poni_name')
        poni = calibration_folder + '\\' + poni_name + '.poni'
        mask_name = self.config.get(section_name, 'mask_name')
        mask = calibration_folder + '\\' + mask_name + '.edf'
        
        return working_folder, poni, mask

    def read_metadata_names(self):
        section_name = 'Image_metadata'
        
        temperature_name_1 = self.config.get(section_name, 'temperature_name_1')
        temperature_name_2 = self.config.get(section_name, 'temperature_name_2')
        temperature_name = [temperature_name_1, temperature_name_2]
        return temperature_name
    
    def read_image_settings(self):
        section_name = 'Image_settings'
        
        colors_list = self.config.get(section_name, 'colors').split(sep=', ')
        colors_list = list([str(name) for name in colors_list])
        color_breakpoints = self.config.get(section_name, 'color_beakpoints').split(sep=', ')
        color_breakpoints = np.array([float(value) for value in color_breakpoints])
        if len(colors_list)!=len(color_breakpoints):
            raise Exception(
                'Error reading the configuration file:\nNumber of colors for '
                +'the colormap ('+str(len(colors_list))+') does not match the '
                +'number of breakpoints ('+str(len(color_breakpoints))+')')
        else:
            color_breakpoint = np.interp(
                color_breakpoints, (color_breakpoints.min(), color_breakpoints.max()),
                (0., 1.))
            custom_cmap = colors.LinearSegmentedColormap.from_list(
                "my_cmap", list(zip(color_breakpoint, colors_list)))
            
        unit = self.config.get(section_name, 'unit')
        npt = int(self.config.get(section_name, 'number_of_points'))
        SUM_IMAGES = eval(self.config.get(section_name, 'SUM_IMAGES'))
        return custom_cmap, unit, npt, SUM_IMAGES

    def get_range_tuple(self, section_name, range_name):
        range_string_list = self.config.get(section_name, range_name).split(sep=', ')
        range_tuple = tuple([float(value) for value in range_string_list])
        return range_tuple

    def read_integ_ranges(self):
        section_name = 'Azimuth_integration_ranges'
        
        azimuth_range_oop = self.get_range_tuple(section_name, 'azimuth_range_oop')
        azimuth_range_ip = self.get_range_tuple(section_name, 'azimuth_range_ip')
        azimuth_range_full = self.get_range_tuple(section_name, 'azimuth_range_full')
        return azimuth_range_oop, azimuth_range_ip, azimuth_range_full

    def read_plotting_ranges(self):
        section_name = 'Plotting_ranges'
        
        int1D_qz_range = self.get_range_tuple(section_name, 'int1D_qz_range')
        int1D_qxy_range = self.get_range_tuple(section_name, 'int1D_qxy_range')
        map2D_qz_range = self.get_range_tuple(section_name, 'map2D_qz_range')
        map2D_qxy_range = self.get_range_tuple(section_name, 'map2D_qxy_range')
        return int1D_qz_range, int1D_qxy_range, map2D_qz_range, map2D_qxy_range


class ConversionFunctions:
    """Big class that contains all the functions needed for the conversion of
    the raw data to the final q_images."""
    
    def direct_beam(self, distance, poni1, poni2, rot1, rot2):
        """Computes the direct beam position (beam center) in pixel coordinates
        from PONI, rotations, and distance."""
        PixelSize=88.54*10**-6 #m
        beamy = (poni1 + np.tan(rot2)/np.cos(rot1) * distance) / PixelSize
        beamx = (poni2 - np.tan(rot1) * distance) / PixelSize
        return beamx,beamy
    
    def SDD(self, distance, rot1, rot2):
        """Returns the effective sample-to-detector distance corrected for
        detector tilts."""
        corrected_distance = distance / (np.cos(rot2) * np.cos(rot1))
        return corrected_distance
    
    def map_pixels(self, shape):
        """Builds 2D pixel coordinate grids for a given image shape."""
        pixelsx = shape[1]
        pixelsy = shape[0]
        x = np.linspace(0, pixelsx - 1, pixelsx)
        y = np.linspace(0, pixelsy - 1, pixelsy)
        map_x, map_y = np.meshgrid(x, y)
        return map_x, map_y
    
    def map_angles(self, map_x, map_y, distance, pixelx, pixely, alpha_i, rot1, rot2, beamx, beamy):
        """Converts pixel positions into scattering angles (theta, alpha)
        and 2*theta using geometry and rotations."""
        X = map_x - beamx
        Y = beamy - map_y
        X *= pixelx #X matrix in mm
        Y *= pixely #Y matrix in mm
        
        #Tilts corrections
        distance = distance / (np.cos(rot2) * np.cos(rot1)) 
        Xp = X * np.cos(rot1)
        Xpp = X * np.sin(rot1)
        Yp = Y * np.cos(rot2)
        Ypp = Y * np.sin(rot2)
        distance = distance - Xpp - Ypp #Distance matrix 
        
        #Conersion in polar units
        r = np.sqrt(Xp ** 2 + Yp ** 2)
        tth = np.arctan2(r, distance)
        phi = np.arctan2(Yp, Xp)
        
        #Conersion in cartesian units
        map_alpha = tth * np.sin(phi) * 180/np.pi - alpha_i
        map_theta = tth * np.cos(phi) * 180/np.pi
        
        return map_theta, map_alpha, tth
    
    def map_q(self, wavelength, map_theta, map_alpha, tth, alpha_i):
        """Converts angles into reciprocal-space components (qx, qy, qz, qr)."""
        map_qx = (2*np.pi/wavelength) * (np.cos(map_theta*np.pi/180) * (
            np.cos(map_alpha*np.pi/180)) - np.cos(alpha_i*np.pi/180))
        map_qy = (2*np.pi/wavelength) * (np.sin(map_theta*np.pi/180) * (
            np.cos(map_alpha*np.pi/180)))
        map_qz = (2*np.pi/wavelength) * (np.sin(map_alpha*np.pi/180)
                                         + np.sin(alpha_i*np.pi/180))
        sign = np.sign(map_qy)
        map_qr = np.sqrt(map_qx**2 + map_qy**2) * sign
        q_par = 4*np.pi/wavelength * np.sin(np.deg2rad(map_theta/2))
        q_per = 4*np.pi/wavelength * np.sin(np.deg2rad(map_alpha/2))
        q = 4*np.pi/wavelength * np.sin(tth/2.)
        return map_qx, map_qy, map_qz, map_qr, q_par, q_per, q

        
class PlottingFunctions:

    def plot_Int1D(self, q, intensity, integ_region, q_limits, filename,
                   sample_metadata, output_root=None):
        """
        Plot and save a 1D intensity curve with metadata in the title.
        
        Parameters
        ----------
        q : 1D numpy.ndarray
            1D array of the q values.
        intensity : 1D numpy.ndarray
            1D array of the integrated scattering intensities.
        integ_region : string
            It classifies the integration region. It should be either OOP, IP
            or ALL.
        q_limits : tuple
            User-defined lower and upper limits for the q axis. 
        filename : string
            Path name of the saved image
        sample_metadata : PlottingMetadata
            Cluster of attributes containng the metadata that will be in the
            final image.
        output_root : string (optional, default None)
            Path of the parent folder of the sample data.
        
        """
        plt.figure(figsize=(8, 5))
        plt.rcParams.update({'font.size': 18})
        plt.plot(q, intensity, 'b-')
        
        plt.title(os.path.basename(sample_metadata.sample_name).replace(
            'rayonix_','')+'_'+str(integ_region)+'\nSPITCH='
            +str(sample_metadata.spitch)+r'$\degree$, T='
            +str(sample_metadata.temperature)+r' $\degree$C, sx= '
            +str(sample_metadata.sx), x=0.5, y=1.05)
        plt.xlim(q_limits[0], q_limits[1])
        plt.ylim(0, 1.1*max(intensity))
        xticks = [q_limits[0]]
        current_tick = 0.5
        while current_tick <= q_limits[1] + 1e-9:
            if current_tick > q_limits[0] + 1e-9:
                xticks.append(current_tick)
            current_tick += 0.5
        plt.xticks(xticks)
        plt.ylabel('I' '(u.a)')
        if integ_region == 'OOP':
            q_label = r'$q_z$($\AA^{-1}$)'
        elif integ_region == 'IP':
            q_label = r'$q_{xy}$($\AA^{-1}$)'
        else:
            q_label = r'$q$($\AA^{-1}$)'
        plt.xlabel(q_label)
        
        output = filename.replace('.edf','_Int1D_'+str(integ_region)+'.jpg')
        if output_root:
            if output_root not in output:
                output = output.replace("RAW", output_root).replace("PROCESSED", output_root)
        else:
            output = output.replace("RAW", "PROCESSED")
        plt.savefig(output, format='jpg', bbox_inches='tight')
        plt.close()
        
    

    def _cell_edges(self, z):
        """Use explicit cell edges to avoid pcolormesh warnings for non-monotonic grids."""
        z_pad = np.pad(z, ((1, 1), (1, 1)), mode='edge')
        return 0.25 * (z_pad[:-1, :-1] + z_pad[1:, :-1] + z_pad[:-1, 1:]
                       + z_pad[1:, 1:])

    def get_color_limits(self, image):
        """Compute suitable min/max intensity limits for plotting on a log scale."""
        cmin = np.round(np.nanmean(image) - 0.1 * np.nanstd(image), 2)
        cmax = np.round(np.nanmean(np.nanmean(image)) + 30* np.nanstd(image), 2)
    
        if cmin < 0:
            cmin = 0.001
        return cmin, cmax
        
    def plot_recmap(self, map_qxy, map_qz, intensity_image, custom_cmap,
                    filename, qxy_limits, qz_limits, sample_metadata, output_root=None,
                    integration_line=None):
        """
        Plot and save a 2D reciprocal-space (qz–qxy) intensity map.
        
        Parameters
        ----------
        map_qxy : 2D numpy.ndarray
            2D matrix of the value of q_xy (= q_r) in nm^-1 for every pixel of
            the detector.
        map_qz : 2D numpy.ndarray
            2D matrix of the value of qz in nm^-1 for every pixel of
            the detector.
        intensity_image : 2D numpy.ndarray
            2D matrix of the rebinned values of normalized scattering intensity.
        custom_cmap : colormap
            User-defined colormap, should fit the experimental peaks of interest.
        filename : string
            Path name of the saved image
        qxy_limits : tuple
            User-defined lower and upper limits for the q_xy axis. 
        qz_limits : tuple
            User-defined lower and upper limits for the q_z axis.
        sample_metadata : PlottingMetadata
            Cluster of attributes containing the metadata that will be in the
            final image.
        output_root : string (optional, default None)
            Path of the parent folder of the sample data.
        integration_line : 1D numpy.ndarray (optional, default None)
            Line delimiting the bounds of integration.
        """
        plt.rcParams.update({'font.size': 20})
        map_qr_edges = self._cell_edges(map_qxy)
        map_qz_edges = self._cell_edges(map_qz)
        color_limits = self.get_color_limits(intensity_image)
        
        plt.pcolormesh(map_qr_edges, map_qz_edges, intensity_image,
                       cmap=custom_cmap, #or turbo, gnuplot,
                       norm=colors.LogNorm(vmin=color_limits[0], vmax=color_limits[1]),
                       shading='auto')
        if integration_line is not None:
            plt.plot_contour(map_qr_edges, map_qz_edges, integration_line)
        
        plt.title(os.path.basename(sample_metadata.sample_name).replace(
            'rayonix_','')+'\nSPITCH='+str(sample_metadata.spitch)
            +r'$\degree$, T='+str(sample_metadata.temperature)
            +r' $\degree$C, sx= '+ str(sample_metadata.sx), x=0.5, y=1.15)
        plt.xlabel(r'$q_{xy}$($\AA^{-1}$)')
        plt.ylabel(r'$q_z$($\AA^{-1}$)')
        
        # Axis limits control the visible q-range in the final qz–qxy plot. (alba)
        plt.xlim(qxy_limits[0], qxy_limits[1])
        #plt.xticks(ticks= (0.0, 0.2, 0.4,0.6), labels=('0.0','0.2','0.4','0.6'))
        plt.ylim(qz_limits[0], qz_limits[1])
        plt.colorbar()
        plt.gca().set_aspect('equal') #Set the aspect ration of the axis scaling
        #plt.tight_layout()
        
        if output_root:
            if output_root not in filename:
               filename = filename.replace("RAW", output_root).replace("PROCESSED", output_root)
        plt.savefig(filename, format='jpg', bbox_inches='tight')
        plt.close()

    
    def plot_azimuthal_map(self, map_tth, map_chi, intensity_image, custom_cmap,
                           filename, sample_metadata, azimuth_range=None,
                           output_root=None):
        """
        Plot and save a 2D reciprocal-space (qz–qxy) intensity map.
        
        Parameters
        ----------
        map_tth : 2D numpy.ndarray
            ???
        map_chi : 2D numpy.ndarray
            2D matrix of the value of the azimuthal angle Chi in deg for every
            pixel of the detector.
        intensity_image : 2D numpy.ndarray
            2D matrix of the rebinned values of normalized scattering intensity.
        custom_cmap : colormap
            User-defined colormap, should fit the experimental peaks of interest.
        filename : string
            Path name of the saved image
        sample_metadata : PlottingMetadata
            Cluster of attributes containing the metadata that will be in the
            final image.
        azimuth_range : tuple (optional, default None)
            Lower and upper limits for the chi axis
        output_root : string (optional, default None)
            Path of the parent folder of the sample data.
        """
          
        color_limits = self.get_color_limits(intensity_image)
        plt.pcolormesh(map_tth, map_chi, intensity_image,
                       cmap=custom_cmap, #turbo, gnuplot,
                       norm=colors.LogNorm(vmin=color_limits[0], vmax=color_limits[1]),
                       shading='auto')
        
        plt.title(os.path.basename(sample_metadata.sample_name).replace(
            'rayonix_','')+'\nSPITCH='+str(sample_metadata.spitch)
            +r'$\degree$, T='+str(sample_metadata.temperature)
            +r' $\degree$C, sx= '+ str(sample_metadata.sx), x=0.5, y=1.15)
        plt.xlabel(r'$q_{xy}$($\AA^{-1}$)')
        plt.ylabel(u'\u03C7'+' ($^o$)')
        
        #plt.xticks(ticks= (0.0, 0.2, 0.4,0.6), labels=('0.0','0.2','0.4','0.6'))
        plt.ylim(0, -180)
        plt.colorbar()
        #plt.gca().set_aspect('equal') #Set the aspect of the axis scaling, i.e. the ratio of y-unit to x-unit.
        #plt.tight_layout()
        #output = File.replace('.edf','_recmap'+'_%.2f_%idegC.jpg'%(sample_metadata.spitch, sample_metadata.temperature))
        if output_root:
            if output_root not in filename:
               filename = filename.replace("RAW", output_root).replace("PROCESSED", output_root)
            
        plt.savefig(filename, format='png', bbox_inches='tight')
        plt.close()
    
    
    def mkdir_filename(self, output):
        """Given an output path, create the corresponding folder."""
        if not os.path.exists(os.path.dirname(output)):
            os.makedirs(os.path.dirname(output))
    
    ##################################################################
    #####                           INFO                         #####
    ##################################################################
    
    def time_info(self, msg):
        """Print the timestamp when called."""
        now = datetime.datetime.now()
        now = now.strftime("%Y-%m-%d %H:%M:%S")
        print("%s --> %s" % (now, str(msg)))

