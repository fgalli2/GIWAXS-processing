##################################################################
#####                                                        #####
#####         IMAGE TREATING PROGRAM FOR GIWAX               #####
#####                                                        #####
#####      Created by Eduardo Solano and Alba Cazorla        #####
#####      Comented by Cristina López-Cava                   #####
#####                                                        #####
##################################################################


import os
import sys
import time #Sleeping
import math, copy

#import cv2
import numpy as np
from PIL import Image #Save tiff file from image
import fabio  # Read/write detector images (e.g., EDF).
import pyFAI  # GIWAXS integration (2D image -> 1D/2D in q-space).
 

from Functions_pyFAI import ReadingInputFunctions, ConversionFunctions, PlottingFunctions
ReadingInput = ReadingInputFunctions()
Conversion = ConversionFunctions()
Plotting = PlottingFunctions()


##############################################################################
class PlottingMetadata:
    def __init__(self, sample_name, temperature, spitch, sx):
        self.sample_name = sample_name
        self.temperature = temperature
        self.spitch = spitch
        self.sx = sx

        
class MakeProfiles():
    def extract_info(self, img):
        filename = img.filename
        sample_name = filename.replace('.edf', '') 
        spitch = float(img.header.get("spitch")) 
        sx = float(img.header.get("sx"))
        
        #To track the fluctuations of the beam
        #MonitorRing = Exposure multiplied by Beam Current
        normalization_factor = float(img.header.get("MonitorRing")) 
        
        if temperature_name[0] in img.header.keys():
            temp_name = temperature_name[0]
        else:
            temp_name = temperature_name[1]
        temperature = round(float(img.header.get(temp_name)),1)
        
        plotting_metadata = PlottingMetadata(sample_name, temperature, spitch, sx)
        
        # Calibration parameters
        mask = self._mask #mask = None
        ai = self._ai
        poni1, poni2, rot1, rot2, _ = ai.poni1, ai.poni2, ai.rot1, ai.rot2, ai.rot3
        wavelength = ai.wavelength * 10**10
        distance_metadata = ai.get_dist()
        beamx, beamy = Conversion.direct_beam(distance_metadata, poni1, poni2,
                                              rot1, rot2)
        distance = Conversion.SDD(distance_metadata, rot1, rot2)
        
        shape = img.shape
        pixelx = pixely = 88.54*10**-6 # in m
                  
    
        
        ##################################################################
        #####                  DATA CONVERSION                      ######
        ##################################################################
        # 1.Pixel map
        map_x, map_y = Conversion.map_pixels(shape)
        
        # 2. Angle map                         
        map_theta, map_alpha,tth = Conversion.map_angles(
            map_x, map_y, distance, pixelx, pixely, spitch, rot1, rot2, beamx,
            beamy)
        
        # 3. q map
        map_q = Conversion.map_q(wavelength, map_theta, map_alpha, tth, spitch)
        map_qz, map_qxy, q = map_q[2], map_q[3], map_q[6]

        
        ##################################################################
        #####                    INTEGRATION                        ######
        ##################################################################
        
        dirname = os.path.dirname(img.filename)
        basename = os.path.basename(img.filename)
             

        ## INTEGRATION ALONG ALL AZIMUTHAL RANGE
        output = os.path.join(dirname, "Int1D_all", basename).replace(
            '.edf','_Int1D_all_%.2f_%idegC.dat' %(spitch, temperature)).replace(
                "RAW", "PROCESSED")
        q_all, I_all = ai.integrate1d(
            data=img.data,
            npt=npt,
            filename=output,
            azimuth_range=(-180,0),
            mask=mask,
            polarization_factor=0.99, 
            method='splitpixel',
            unit=unit,
            normalization_factor=normalization_factor)
    
        Plotting.mkdir_filename(output)
        output= os.path.join(dirname, "Int1D_all", basename).replace(
            '.edf','_Int1D_'+'all'+'_%.2f_%idegC.jpg' %(spitch, temperature)).replace(
                "RAW", "PROCESSED")
        Plotting.plot_Int1D(
            q=q_all,
            intensity=I_all,
            integ_region='ALL',
            q_limits=(np.min(q_all),np.max(q_all)),
            filename=output,
            sample_metadata=plotting_metadata)
            
  
        ## OUT OF PLANE INTEGRATION 
        
        output = os.path.join(dirname, "oop", basename).replace(
            ".edf", '_q_A^-1_Int_OOP_%.2f_%idegC.dat' %(spitch, temperature)).replace(
                "RAW", "PROCESSED")
        q, I_OOP = ai.integrate1d(
            data=img.data,
            npt=npt,
            filename=output,
            azimuth_range=azimuth_range_oop,
            mask=mask,
            polarization_factor=0.99,
            method='splitpixel', 
            unit=unit,
            normalization_factor=normalization_factor)
        
        output = os.path.join(dirname, "oop", basename).replace(
            '.edf','_Int1D_'+'OOP'+'_%.2f_%idegC.jpg' %(spitch, temperature)).replace(
                "RAW", "PROCESSED")
        Plotting.mkdir_filename(output)
        Plotting.plot_Int1D(
            q=q,
            intensity=I_OOP,
            integ_region='OOP',
            q_limits=int1D_qz_range,
            filename=output,
            sample_metadata=plotting_metadata)
        
     
        ## INTEGRATION ALONG THE IN-PLANE DIRECTION 
        
        output = os.path.join(dirname, "ip", basename).replace(
            ".edf", '_q_A^-1_ip_%.2f_%idegC.dat' %(spitch, temperature)).replace(
                "RAW", "PROCESSED")
        q, I_IP = ai.integrate1d(
            data=img.data,
            npt=npt,
            filename=output,
            azimuth_range=azimuth_range_ip,
            mask=mask,
            polarization_factor=0.99,
            method='splitpixel', 
            unit=unit,
            normalization_factor=normalization_factor)
        
        output = os.path.join(dirname, "ip", basename).replace(
            '.edf','_Int1D_'+'IP'+'_%.2f_%idegC.jpg' %(spitch, temperature)).replace(
                "RAW", "PROCESSED")
        Plotting.mkdir_filename(output)
        Plotting.plot_Int1D(
            q=q,
            intensity=I_IP,
            integ_region='IP',
            q_limits=int1D_qxy_range,
            filename=output,
            sample_metadata=plotting_metadata)
        
        
        ## 2D GIWAXS MAPS (qz–qxy reciprocal-space maps) ##
        intensity_image = copy.copy(img.data)
        
        #nan pixels along qz for qxy approx 0 since we are not probing
        #the true specular direction
        intensity_image[:, math.floor(beamx)-2:math.ceil(beamx) + 1] = np.nan 
        output = os.path.join(dirname, "recmap", basename).replace(
            '.edf','_recmap'+'_%.2f_%idegC.jpg'%(spitch, temperature)).replace(
                '.edf', '_recmap.jpg').replace("RAW", "PROCESSED")
        Plotting.mkdir_filename(output)
        Plotting.plot_recmap(
            map_qxy=map_qxy,
            map_qz=map_qz,
            intensity_image=intensity_image,
            custom_cmap=custom_cmap,
            filename=output,
            qxy_limits=map2D_qxy_range,
            qz_limits=map2D_qz_range,
            sample_metadata=plotting_metadata)


        ## 2D INTEGRATION OVER AZMIIMUTH        
        output = os.path.join(dirname, "Int2D_azimuth", basename).replace(
            '.edf','_Int2Dazimuth'+'_%.2f_%idegC.tiff' %(spitch, temperature)).replace(
                "RAW", "PROCESSED")
        
        intensity_image_az, map_tth, map_chi = ai.integrate2d(
            data=img.data,
            filename=output,
            azimuth_range=(-180,0),
            npt_rad=1000,
            npt_azim=1000,
            polarization_factor=0.99,
            method='splitpixel',
            unit=unit,
            normalization_factor=normalization_factor) 
        intensity_image_az_tiff = Image.fromarray(
            intensity_image_az.astype(np.float32), mode='F')
        intensity_image_az_tiff.save(output, 'TIFF')
        
        Plotting.mkdir_filename(output)
        output = os.path.join(dirname, "Int2D_azimuth", basename).replace(
            '.edf','_Int2Dazimuth'+'_%.2f_%idegC.png' %(spitch, temperature)).replace(
                "RAW", "PROCESSED")
        Plotting.plot_azimuthal_map(
            map_tth,
            map_chi,
            intensity_image_az,
            custom_cmap,
            output,
            sample_metadata=plotting_metadata)


##################################################################
#####                     IMAGE SUMMATION                   #####
##################################################################

class SumImages():
    ADD_CHANNELS = ['Monitor', 'MonitorRing', 'Photo', 'exposure',
                    'acq_expo_time', 'tfg_ch1', 'tfg_ch2', 'tfg_ch3',
                    'tfg_ch4', 'tfg_ch5', 'tfg_ch6', 'tfg_ch7', 'tfg_ch8']
    EXT = ".edf"
    OUTPUT_NAME = "_add"

    def sum_images(self, filename):
        base_edf = filename[:-8]
        basename = os.path.basename(base_edf)
        dirname = os.path.dirname(base_edf)
        image_group = [os.path.join(dirname, f) for f in os.listdir(dirname)
                       if basename in f and filename not in f]

        if len(image_group) < 2:
            return

        base_image = fabio.open(filename) #The first image is taken as the base
        base_image.data = base_image.data.astype(np.float32)
        for image in image_group:
            Plotting.time_info('Combining image "%s"' % image)
            new_image = fabio.open(image)
            new_image.data = new_image.data.astype(np.float32)
            for ch in self.ADD_CHANNELS:
                try:
                    base_image.header[ch] = float(
                        new_image.header[ch]) + float(base_image.header[ch])
                except:
                    pass
            base_image.data += new_image.data

        output = filename.replace("RAW", "PROCESSED").replace(
            self.EXT, self.OUTPUT_NAME+self.EXT)
        if not os.path.exists(os.path.dirname(output)):
            os.makedirs(os.path.dirname(output))
        base_image.save(output)
        Plotting.time_info('Images combined and saved as: "%s"' % output)
        return output

    
    
##################################################################
#####                       MAIN PROCESS                    #####
##################################################################

class Process(MakeProfiles):
    def __init__(self):
        Plotting.time_info("Preparing environment")
        self.ai = pyFAI.load(poni)
        self.mask = fabio.open(mask).data
        self._files_integrated = list()

    def get_files2integrate(self, folder):
        Plotting.time_info("Processing folder: %s" % folder)
        files2integrate = [f for f in self.get_files(folder=folder)
                           if f not in self._files_integrated]
        Plotting.time_info("%s file(s) found to be integrated"
                           % len(files2integrate))
        return files2integrate

    def get_files(self, folder):
        files = [os.path.join(folder, f) for f in os.listdir(folder)
                 if f.endswith(".edf")]
        Plotting.time_info('Initial files:\n  '+'\n  '.join(files))
        if SUM_IMAGES:
            files = [f for f in files if "0000.edf" in f]
            Plotting.time_info('Files after addition:\n  '+'\n  '.join(files))
        return files

    def process_folder(self, folder):
        files2integrate = self.get_files2integrate(folder=folder)
        if not files2integrate: return
        time.sleep(1)

        for i, filename in enumerate(files2integrate):
            filename_base = filename
            if SUM_IMAGES:
                filename_base = SumImages().sum_images(filename)
                if filename_base is None:
                    self._files_integrated.append(filename)
                    continue
            self.process_filename(filename_base)
            self._files_integrated.append(filename)

    
    def select_ai(self, edf):
        self._ai = self.ai
        self._mask = self.mask
        
    def ignore_edf(self, edf):
        return False


    def process_filename(self, filename):
        edf = fabio.open(filename)
        
        if self.ignore_edf(edf):
            return
        self.select_ai(edf)
        self.extract_info(edf)
        self._files_integrated.append(filename)
        print("File successfully processed: %s" % edf.filename)



##################################################################
#####                         ENTRYPOINT                    #####
##################################################################

if __name__=="__main__":
    folder, poni, mask = ReadingInput.read_input_files()
    temperature_name = ReadingInput.read_metadata_names()
    custom_cmap, unit, npt, SUM_IMAGES = ReadingInput.read_image_settings()
    (azimuth_range_oop, azimuth_range_ip,
     azimuth_range_full) = ReadingInput.read_integ_ranges()
    (int1D_qz_range, int1D_qxy_range, map2D_qz_range,
     map2D_qxy_range) = ReadingInput.read_plotting_ranges()
    
    Process().process_folder(folder)
    sys.exit()

