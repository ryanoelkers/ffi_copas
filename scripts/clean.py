""" This script will clean the FFI images, it has been significantly reduced to only
show specific steps. The majority of the functional library can be found in the libraries
directory in utils, fits, photometry, image, and preprocessing libraries."""
from libraries.utils import Utils
from libraries.preprocessing import Preprocessing
from scripts.master import Master
from config import Configuration
import os
from astropy.io import fits
import time
import numpy as np


class Clean:

    @staticmethod
    def clean_images(sky_subtract='Y', bias_subtract="N", flat_divide='N', alignment='N'):
        """ This is the main function script to clean multiple images, alternatively clean_img can be used to clean
        a single image.
        :parameter sky_subtract - Y/N if you want to subtract the sky from the images (default = Y)
        :parameter bias_subtract - Y/N if you want to subtract the bias from the images (default = N)
        :parameter flat_divide - Y/N if you want to flatten the images (default = N)
        :parameter alignment - Y/N if you want to align the images (default = N)

        return no value is returned, the values images from in_path are cleaned and deposited in out_path
        """
        st = time.time()  # clock started

        # get the file list
        Utils.log("Getting file list...", "info", Configuration.LOG_SCREEN)
        files = Utils.get_file_list(Configuration.RAW_DIRECTORY, '-' + Configuration.CAMERA +
                                    '-' + Configuration.CCD + '-' + Configuration.SECT_NUM + '-s_ffic.fits')

        # loop through the images, look for the OK images, and get the image to be used for the master frame
        img_chk, ref_path = Master.check_frame_quality(files)

        # break if there are no files
        if len(files) == 0:
            Utils.log("No .fits files found in " + Configuration.RAW_DIRECTORY + ". Breaking...",
                      "debug", Configuration.LOG_SCREEN)
            return()
        
        Utils.log("Starting to clean " + str(len(files)) + " images.", "info", Configuration.LOG_SCREEN)
        for idx, file in enumerate(files):

            # make a new name for the file based on which actions are taken
            file_name = Preprocessing.mk_nme(file, 'N', sky_subtract, bias_subtract, flat_divide, alignment)

            # only create the files that don't exist
            if os.path.isfile(Configuration.CLEAN_DIRECTORY + file_name) == 1:
                Utils.log("Image " + Configuration.CLEAN_DIRECTORY + file_name + " already exists. Skipping for now...",
                          "info", Configuration.LOG_SCREEN)

            # if the image does not exist then clean
            if os.path.isfile(Configuration.CLEAN_DIRECTORY + file_name) == 0:

                # clean the image
                clean_img, header, bd_flag = Clean.clean_img(file, ref_path, sky_subtract,
                                                             bias_subtract, flat_divide, alignment)

                # write out the file
                if bd_flag == 0:
                    fits.writeto(Configuration.CLEAN_DIRECTORY + file_name, clean_img, header, overwrite=True)

                    # print an update to the cleaning process
                    Utils.log("Cleaned image written as " + Configuration.CLEAN_DIRECTORY + file_name + ".",
                              "info", Configuration.LOG_SCREEN)
                else:
                    Utils.log(file_name + " is a bad image. Not written.", "info", Configuration.LOG_SCREEN)

            Utils.log(str(len(files) - idx - 1) + " images remain to be cleaned.",  "info", Configuration.LOG_SCREEN)

        fn = time.time()  # clock stopped
        Utils.log("Imaging cleaning complete in " + str(np.around((fn - st), decimals=2)) + "s.", "info", 'Y')

    @staticmethod
    def clean_img(file, ref_path, sky_subtract='Y', bias_subtract='N', flat_divide='N', alignment='N'):
        """ This function is the primary script to clean the image, various other functions found in this class
        can be found in the various libraries imported.

        :parameter  file - The file name of the image you would like to clean
        :parameter ref_path - The path to the reference frame
        :parameter sky_subtract - Y/N if you want to subtract the sky from the image (default = Y)
        :parameter bias_subtract - Y/N if you want to remove a bias frame (default = N)
        :parameter flat_divide - Y/N if you want to flatten the image (default = N)
        :parameter alignment - Y/N if you want to align the image (default = N)
        """

        Utils.log("Now cleaning " + file + ".", "info", 'Y')

        # read in the image
        org_img, org_header = fits.getdata(Configuration.RAW_DIRECTORY + file, header=True)

        if np.median(org_img) > 0:
            # remove the overscan regions
            big_img, header = Preprocessing.remove_overscan(org_img, org_header)

            # bias subtract if necessary
            if bias_subtract == 'Y':
                st = time.time()
                big_img, header = Preprocessing.bias_subtract(big_img, header)
                fn = time.time()
                Utils.log("Image bias corrected in " + str(np.around((fn - st), decimals=2)) + "s.",
                          "info", Configuration.LOG_SCREEN)

            if bias_subtract == 'N':
                Utils.log("Skipping bias correction....", "info", Configuration.LOG_SCREEN)

            # flat divide if necessary
            if flat_divide == 'Y':
                st = time.time()
                big_img, header = Preprocessing.flat_divide(big_img, header)
                fn = time.time()
                Utils.log("Image flattened in " + str(np.around((fn - st), decimals=2)) + "s.",
                          "info", Configuration.LOG_SCREEN)

            if flat_divide == 'N':
                Utils.log("Skipping image flattening....", "info", Configuration.LOG_SCREEN)

            # sky subtract if necessary
            if sky_subtract == 'Y':
                st = time.time()
                # the background sample size is set to pix x pix pixels, and bxs x bxs sub images
                # this should not be hard coded...update for later

                Utils.log("A background box of " + str(Configuration.PIX) + " x " + str(Configuration.PIX) +
                          " will be used for background subtraction.", "info", Configuration.LOG_SCREEN)

                big_img, header = Preprocessing.sky_subtract(big_img, header, 'N')
                fn = time.time()
                Utils.log("Sky subtracted in " + str(np.around((fn - st), decimals=2)) + "s.",
                          "info", Configuration.LOG_SCREEN)

            if sky_subtract == 'N':
                Utils.log("Skipping sky subtraction...", "info", Configuration.LOG_SCREEN)

            # align the image if necessary
            if alignment == 'Y':
                st = time.time()
                big_img, header = Preprocessing.align_img(big_img, header, ref_path)
                fn = time.time()
                Utils.log("Image aligned in " + str(np.around((fn - st), decimals=2)) + "s.",
                          "info", Configuration.LOG_SCREEN)

            if alignment == 'N':
                Utils.log("Skipping image alignment....", "info", Configuration.LOG_SCREEN)

            Utils.log("Cleaning finished.", "info", Configuration.LOG_SCREEN)
            bd_flag = 0
        else:
            Utils.log("Bad image!", "info", Configuration.LOG_SCREEN)
            big_img = org_img
            header = org_header
            bd_flag = 1

        return big_img, header, bd_flag
