""" This script will generate the master frame (either the first frame in a series, or a combinatoin of frames),
as well as generate the star lists for the field"""
from libraries.utils import Utils
from libraries.dbaccess import DBaccess
from libraries.images import Images
from libraries.photometry import Photometry
from config import Configuration
import os
import pandas as pd
import numpy as np
from astropy.io import fits


class Master:

    @staticmethod
    def multi_master(img_chk):
        """ This function will generate a master frame from multiple images. It uses a median of the available images,
        and only combines those with quality flags of 0. This function works in bulks of 100 images for combination.

        :parameter img_chk - The data frame with the image information

        :return master_image - The master image which is based on the combination
        :return master_header - The header of the master image, only the WCS coordinates are kept in the image
        """

        if os.path.isfile(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" + Configuration.CAMERA + "_" +
                          Configuration.CCD + "_master.fits") == 0:
            # remove all of the images from the image check data frame with bad quality flags
            img_good = img_chk[img_chk['pass'] == 0].copy().reset_index(drop=True)

            # determine the number of loops we need to move through for each image
            nfiles = len(img_good)
            nbulk = 20

            # get the integer and remainder for the combination
            full_bulk = nfiles // nbulk
            part_bulk = nfiles % nbulk
            if part_bulk > 0:
                hold_bulk = full_bulk + 1
            else:
                hold_bulk = full_bulk

            # here is the 'holder'
            hold_data = np.ndarray(shape=(hold_bulk, Configuration.AXS, Configuration.AXS))

            # update the log
            Utils.log("Generating a master frame from multiple files in bulks of " + str(nbulk) +
                      " images. There are " + str(nfiles) +
                      " good images to combine, which means there should be " + str(hold_bulk) +
                      " mini-files to combine.", "info", "Y")

            for kk in range(0, hold_bulk):

                # loop through the images in sets of nbulk
                if kk < full_bulk:
                    # generate the image holder
                    block_hold = np.ndarray(shape=(nbulk, Configuration.AXS, Configuration.AXS))

                    # generate the max index
                    mx_index = nbulk
                else:
                    # generate the image holder
                    block_hold = np.ndarray(shape=(part_bulk, Configuration.AXS, Configuration.AXS))

                    # generate the max index
                    mx_index = part_bulk

                # make the starting index
                loop_start = kk * nbulk
                idx_cnt = 0

                Utils.log("Making mini-file " + str(kk) + ".....", "info", "Y")

                # now loop through the images
                for jj in range(loop_start, mx_index + loop_start):
                    # read in the image directly into the block_hold
                    block_hold[idx_cnt] = fits.getdata(Configuration.CLEAN_DIRECTORY + img_good['file'][jj])

                    # increase the iteration
                    idx_cnt += 1

                # median the data into a single file
                hold_data[kk] = np.median(block_hold, axis=0)

            # median the mini-images into one large image
            master_image = np.median(hold_data, axis=0)

            # pull the header information from the first file of the set
            master_header = fits.getheader(Configuration.CLEAN_DIRECTORY +
                                           img_chk[img_chk['pass'] == 0].reset_index(drop=True)['file'][0])
            master_header['master'] = 'multi'
            master_header['num_comb'] = nfiles

            # write the image out to the master directory
            fits.writeto(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" + Configuration.CAMERA + "_" +
                         Configuration.CCD + "_master.fits", master_image, master_header, overwrite=True)

            Utils.log("Multi-master frame created and written to master directory.", "info", "Y")
        else:
            Utils.log("Legacy multi-master frame found. Reading in legacy data.", "info", "Y")

            master_image = fits.getdata(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" +
                                        Configuration.CAMERA + "_" + Configuration.CCD + "_master.fits")

            master_header = fits.getheader(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" +
                                           Configuration.CAMERA + "_" + Configuration.CCD + "_master.fits")

        return master_image, master_header

    @staticmethod
    def check_frame_quality(file_list):
        """ This function will check the individual quality flags to determine which files should be used for the master
        frame. It will print out the data to the master frame directory, in order to save information for later

        :parameter file_list - The list of files to investigate for quality
        :return img_chk - A data frame with the quality information for each frame
        :return ref_path - The path to the first file in the series which will be the reference frame
        """

        if os.path.exists(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" +
                          Configuration.CAMERA + "_" + Configuration.CCD + "_image_quality.csv"):
            Utils.log("Legacy image quality information found. Reading in legacy file.", "info", "Y")

            # read in the legacy quality information
            img_chk = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" +
                                  Configuration.CAMERA + "_" + Configuration.CCD + "_image_quality.csv", index_col=0)
        else:
            Utils.log("Determining image quality information.", "info", "Y")
            # set up a variable to look to 'pass' the image based on data quality
            img_chk = pd.DataFrame(columns=['file', 'JD', 'ra', 'dec', 'pass'])

            # loop through reach file in the file_list and look for the quality flag of 0
            for idx, file in enumerate(file_list):

                if fits.getheader(Configuration.RAW_DIRECTORY + file, 1)['DQUALITY'] == 0:
                    # sometimes the quality flag does not catch the issue
                    try:
                        # add the data to a data frame which will be returned for analysis
                        img_chk = img_chk.append(pd.DataFrame(
                            data={'file': [file],
                                  'JD': [np.mean([fits.getheader(Configuration.RAW_DIRECTORY + file, 1)['TSTART'],
                                                 fits.getheader(Configuration.RAW_DIRECTORY + file, 1)['TSTOP']])],
                                  'ra': [fits.getheader(Configuration.RAW_DIRECTORY + file, 1)['CRVAL1']],
                                  'dec': [fits.getheader(Configuration.RAW_DIRECTORY + file, 1)['CRVAL2']],
                                  'pass': [fits.getheader(Configuration.RAW_DIRECTORY + file, 1)['DQUALITY']]})
                        ).reset_index(drop=True)
                    except KeyError:
                        img_chk = img_chk.append(pd.DataFrame(
                            data={'file': [file],
                                  'JD': [np.mean([fits.getheader(Configuration.RAW_DIRECTORY + file, 1)['TSTART'],
                                                  fits.getheader(Configuration.RAW_DIRECTORY + file, 1)['TSTOP']])],
                                  'ra': -99,
                                  'dec': -99,
                                  'pass': 99})
                        ).reset_index(drop=True)
                else:
                    # add the data to a data frame which will be returned for analysis
                    img_chk = img_chk.append(pd.DataFrame(
                        data={'file': [file],
                              'JD': [np.mean([fits.getheader(Configuration.RAW_DIRECTORY + file, 1)['TSTART'],
                                             fits.getheader(Configuration.RAW_DIRECTORY + file, 1)['TSTOP']])],
                              'ra': -99,
                              'dec': -99,
                              'pass': [fits.getheader(Configuration.RAW_DIRECTORY + file, 1)['DQUALITY']]})
                    ).reset_index(drop=True)

            # dump the file to make a legacy file for later
            img_chk.to_csv(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" +
                           Configuration.CAMERA + "_" + Configuration.CCD + "_image_quality.csv")

        # get the path to the file with the first OK frame. this is for the reference positioning later
        ref_path = Configuration.RAW_DIRECTORY + img_chk[img_chk['pass'] == 0].reset_index(drop=True)['file'][0]

        return img_chk, ref_path

    @staticmethod
    def mk_master(mast_type='single'):
        """ This function will generate the master frame and generates position files, or if a single frame is chosen,
        then only the position file are generated.

        :parameter mast_type - Default is the first frame in the series ('single'), 'combine' will generate from all
                                CCDs in the directory
        :return - Either the master frame is generated, or the single frame is copied over and a star list is generated
        """

        if (os.path.isfile(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" +
                           Configuration.CAMERA + '_' +
                           Configuration.CCD + '_master.data') == 1) & \
                (os.path.isfile(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" +
                                Configuration.CAMERA + "_" + Configuration.CCD + "_ticv7_list.csv") == 1) & \
                (os.path.isfile(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" +
                                Configuration.CAMERA + '_' + Configuration.CCD + '_master.fits') == 1):

            Utils.log("Legacy master frame information found, pulling legacy master frame information.", "info", 'Y')

            # read in the legacy data
            master_df = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" +
                                    Configuration.CAMERA + '_' + Configuration.CCD + '_master.data', index_col=0)

        if (os.path.isfile(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" +
                           Configuration.CAMERA + '_' +
                           Configuration.CCD + '_master.data') == 0) | \
            (os.path.isfile(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" + Configuration.CAMERA + "_" +
                            Configuration.CCD + "_ticv7_list.csv") == 0) | \
                (os.path.isfile(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" +
                                Configuration.CAMERA + '_' + Configuration.CCD + '_master.fits') == 0):

            # get the file list
            Utils.log("Master frame information not found, getting master frame, and photometry files.", "info", 'Y')
            files = Utils.get_file_list(Configuration.CLEAN_DIRECTORY, '-' + Configuration.CAMERA +
                                        '-' + Configuration.CCD + '-' + Configuration.SECT_NUM + '-s_ffic_' +
                                        Configuration.FILE_EXT + '.fits')

            img_chk, ref_path = Master.check_frame_quality(files)
            img_chk['file'] = img_chk.apply(lambda x:
                                            x['file'].split('.')[0] + '_' + Configuration.FILE_EXT + '.' + x['file'].
                                            split('.')[1], axis=1)

            # now separate depending on whether there is a single master or combination master
            if mast_type == 'single':
                Utils.log("Generating a master frame from a single file. Using the first available good quality frame.",
                          "info", "Y")
                # read in the first frame from the directory and then read it in
                master = fits.getdata(Configuration.CLEAN_DIRECTORY +
                                      img_chk[img_chk['pass'] == 0].reset_index(drop=True)['file'][0])
                mast_head = fits.getheader(Configuration.CLEAN_DIRECTORY +
                                           img_chk[img_chk['pass'] == 0].reset_index(drop=True)['file'][0])
                mast_head['master'] = 'single'
                mast_head['num_comb'] = 1

                # write the image out to the master directory
                fits.writeto(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" + Configuration.CAMERA + "_" +
                             Configuration.CCD + "_master.fits", master, mast_head, overwrite=True)
            else:
                # generate the master frame from the images which have a quality flag set to 0
                master, mast_head = Master.multi_master(img_chk)

            # get the center pixels of the master frame, and the maximum angular distance to the edge
            mx_dist, cen_ra, cen_de = Images.get_coords(mast_head)

            if Configuration.MACHINE == 'tessdev':
                # read in the appropriate sql file, then replace the necessary values
                sql_cmd = DBaccess.get_starlist_query(Configuration.QUERIES_DIRECTORY + 'star_list.sql',
                                                      cen_ra, cen_de, mx_dist)

                # now query the TIC at this given distance and center coordinate position
                tic_list = DBaccess.query_tic7(sql_cmd, Configuration.MASTER_DIRECTORY,
                                               Configuration.SECTOR + "_" + Configuration.CAMERA + "_" +
                                               Configuration.CCD + "_ticv7_list.csv", Configuration.MACHINE)
            else:
                tic_list = DBaccess.query_mast(cen_ra, cen_de, mx_dist,
                                               Configuration.MASTER_DIRECTORY,
                                               Configuration.SECTOR + "_" + Configuration.CAMERA + "_" +
                                               Configuration.CCD + "_ticv7_list.csv")
            Utils.log("Running photometry on master frame.", "info", 'Y')

            # now we run photometry on the master frame
            phot_df = Photometry.single_frame_aper(master, tic_list, mast_head, master_frame='Y')

            # remove the unseen stars
            phot_df = phot_df.dropna()

            # clip stars outside of the frame
            master_df = phot_df[(phot_df.x > 0) & (phot_df.x < Configuration.AXS) &
                                (phot_df.y > 0) & (phot_df.y < Configuration.AXS) &
                                (phot_df.mag_err < Configuration.ERROR_LIMIT)].copy()

            # sort the list by magnitude
            master_df = master_df.sort_values(by='tessmag', ascending=1).reset_index(drop=True)
            Utils.log(str(len(master_df)) + " stars found with OK photometry. Dumping master frame files for backup.",
                      "info", 'Y')

            # dump the various photometry files for future reference
            master_df[['TICID', 'tessmag', 'ra', 'dec',
                       'x', 'y', 'mag', 'mag_err', 'flux', 'flux_err']].to_csv(Configuration.MASTER_DIRECTORY +
                                                                               Configuration.SECTOR + "_" +
                                                                               Configuration.CAMERA + '_' +
                                                                               Configuration.CCD + '_master.data')

        return master_df
