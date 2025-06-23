from astropy.io import fits
import pandas as pd
import numpy as np
from astropy.stats import sigma_clipped_stats
from config import Configuration
from libraries.utils import Utils
import os

class Detection:

    @staticmethod
    def mk_mask(radius, pixs):
        """ This function will generate the masks for specific Sector, Camera, CCD combinations. If a mask already
        exists, then a new one is not generated.

        :parameter radius - The radius used in photometry from Oelkers+2018, this is the size of the stellar mask
        :parameter pixs - How many pixels (on a side) is the data?

        :return mask_img - The mask image is returned.
        """

        if os.path.isfile(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + '_' +
                          Configuration.CAMERA + '_' + Configuration.CCD + "_mask.fits"):
            Utils.log("A legacy mask has been found for Sector: " + Configuration.SECTOR + " Camera: " +
                      Configuration.CAMERA + " CCD: " + Configuration.CCD +
                      ". If you do not want to use this mask, delete the file!", "info", Configuration.LOG_SCREEN)
            # pull in the legacy file if it exists
            mask_img = fits.getdata(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + '_' +
                                    Configuration.CAMERA + '_' + Configuration.CCD + "_mask.fits", header=False)

        else:
            Utils.log("No legacy mask found for Sector: " + Configuration.SECTOR + " Camera: " + Configuration.CAMERA +
                      " CCD: " + Configuration.CCD + ". Generating mask now.", "info", Configuration.LOG_SCREEN)

            # make the mask using the aperture file
            mask_img = np.zeros((pixs, pixs))

            # get the star list for the data
            star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + '_' +
                                    Configuration.CAMERA + '_' + Configuration.CCD + '_master.data',
                                    sep=',', header=0, index_col=0)

            # # make a mesh grid for the mask
            y, x = np.mgrid[0:pixs, 0:pixs]

            # # go through each star in the star list
            for idx, star in star_list.iterrows():

                # make the mask position
                mask2 = (radius ** 2 > (x - star.x) ** 2 + (y - star.y) ** 2)
                mask_img[mask2] = 1

                if idx % 1000 == 0:
                    Utils.log("Working on generating the mask for the next 1000 stars. " + str(len(star_list) - idx) +
                              " stars remain.", "info", Configuration.LOG_SCREEN)

            # write out the image mask
            fits.writeto(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + '_' +
                         Configuration.CAMERA + '_' + Configuration.CCD + "_mask.fits",
                         mask_img, overwrite=True)

        return mask_img

    @staticmethod
    def make_detection_frames(bin_scale):
        """ This function acts as a wrapper function to generate a mask frame and then generate the ensemble of
        detection frames for the current sector, ccd, camera combination.

        :parameter bin_scale - The time frame (in number of frames) you want to bin the detection frames. There will be
                            a final "full-image" detection frame generated which is the entire period binned together."

        :return Nothing is returned but the detection frames are output to a directory
        """

        # get the file list of the differenced images
        file_list = Utils.get_file_list(Configuration.DIFFERENCED_DIRECTORY, '.fits')

        # generate or pull the image mask
        mask_img = Detection.mk_mask(Configuration.APER_SIZE, Configuration.AXS)

        # loop through the differenced images and combine when necessary
        bin_idx = 1
        for idx, file in enumerate(file_list):

            try:
                # get the file
                img, header = fits.getdata(Configuration.DIFFERENCED_DIRECTORY + "/" + file, 0, header=True)

                # get the background statistics on the image
                mn_img, md_img, sd_img = sigma_clipped_stats(img, sigma=2)

                # make the image
                img[mask_img == 1] = np.random.normal(loc=mn_img, scale=sd_img, size=len(np.argwhere(mask_img == 1)))
                img = np.abs(img) / np.sqrt(np.abs(img))

                # co-add full detection frame
                if idx == 0:
                    full_image = img
                    st_time_full = header['TSTART']
                else:
                    full_image = img + full_image
                    ed_time_full = header['TSTOP']

                # co-add the bin_scale frames
                if bin_idx == 1:
                    # add the image to the bin image
                    bin_img = img
                    st_time_bin = header['TSTART']
                    bin_idx = bin_idx + 1
                else:

                    # check the date range before binning
                    ed_time_chk = header['TSTOP']

                    if ed_time_chk - st_time_bin < bin_scale:
                        # if we are still within the time frame, bin and then keep moving
                        bin_img = img + bin_img
                        ed_time_bin = header['TSTOP']

                        # increase bin idx
                        bin_idx = bin_idx + 1

                    else:
                        # get the binning information
                        header['TSTART'] = st_time_bin
                        header['TSTOP'] = ed_time_bin
                        header['NBIN'] = bin_idx
                        tme = int(np.around(np.mean([st_time_bin, ed_time_bin]), decimals=0))

                        nme = Configuration.SECTOR + "_" + \
                              Configuration.CAMERA + "_" + \
                              Configuration.CCD + "_" + str(tme) + "_detection_frame.fits"

                        # write the file
                        fits.writeto(Configuration.DETECTION_FRAME_DIRECTORY + nme, bin_img, header, overwrite='True')
                        Utils.log("Detection frame: " + nme + " has been written.",
                                  "info", Configuration.LOG_SCREEN)
                        bin_idx = 1

            except TypeError:
                Utils.log("Bad file...skipping!", "info", "N")

            if idx % 100 == 0:
                Utils.log("Working on next 100 images. " + str(len(file_list) - idx) + " files remain.",
                          "info", Configuration.LOG_SCREEN)

        # get the proper header information
        header['TSTART'] = st_time_full
        header['TSTOP'] = ed_time_full
        header['NBIN'] = idx

        # write the file
        fits.writeto(Configuration.DETECTION_FRAME_DIRECTORY +
                     Configuration.SECTOR + "_" +
                     Configuration.CAMERA + "_" + Configuration.CCD + "_full_detection_frame.fits",
                     full_image, header, overwrite='True')

        Utils.log("Finished writing detection frames! Full detection frame has been generated.",
                  "info",
                  Configuration.LOG_SCREEN)
        return