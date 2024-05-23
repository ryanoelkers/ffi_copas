""" This script will generate the master frame (either the first frame in a series, or a combinatoin of frames),
as well as generate the star lists for the field"""
from libraries.utils import Utils
from libraries.photometry import Photometry
from config import Configuration
import numpy as np
import os
import pandas as pd
from astropy.io import fits
from scripts.master import Master
from scripts.mast_release import MastRelease
import matplotlib.pyplot as plt


class LightCurves:

    @staticmethod
    def diff_img_phot(star_list):
        """ This fucntion will perform the aperture photometry on the differenced images

        :parameter star_list - The star list provided by the TIC and the master frame

        :return Nothing is returned, but the differenced files are output
        """
        # get the file list of the differenced images
        files = Utils.get_file_list(Configuration.DIFFERENCED_DIRECTORY, '-' + Configuration.SECT + '-' +
                                    Configuration.CAMERA + '-' + Configuration.CCD + '-' +
                                    Configuration.SECT_NUM + '-s_ffic_' + Configuration.FILE_EXT + 'd.fits')
        nfiles = len(files)

        Utils.log(str(nfiles) + ' files found in ' + Configuration.DIFFERENCED_DIRECTORY + ' for Sector: ' +
                  str(Configuration.SECTOR) + ' Camera: ' + Configuration.CAMERA + ' CCD: ' + Configuration.CCD,
                  'info', Configuration.LOG_SCREEN)

        # loop through files, generating the flux files
        for idx, file in enumerate(files):

            # generate the file name
            file_prt = file.split('.')
            flx_file_nme = file_prt[0] + '.flux'

            if os.path.exists(Configuration.DIFFERENCED_DIRECTORY + flx_file_nme) == 1:
                Utils.log(flx_file_nme + " file exists. Skipping...", 'info', Configuration.LOG_SCREEN)

            if os.path.exists(Configuration.DIFFERENCED_DIRECTORY + flx_file_nme) == 0:
                # read in the image
                img, header = fits.getdata(Configuration.DIFFERENCED_DIRECTORY + file, header=True)
                master, mheader = fits.getdata(Configuration.MASTER_DIRECTORY +
                                               Configuration.SECTOR + "_" +
                                               Configuration.CAMERA + "_" +
                                               Configuration.CCD + "_master.fits", header=True)

                # get the JD for the image
                jd1 = header['TSTART']
                jd2 = header['TSTOP']
                jd = np.mean([jd1, jd2])

                Utils.log('Generating flux file for ' + file + '. ' + str(nfiles - idx - 1) +
                          ' files left after this image.', 'info', Configuration.LOG_SCREEN)

                # get photometry for those stars
                star_phot, target_list = Photometry.single_frame_aper(img, star_list, mheader,
                                                                      stars_to_phot=Configuration.STAR_CHECK,
                                                                      bkg_sub='local',
                                                                      flux_only='Y',
                                                                      diff_flux_convert='Y',
                                                                      offset='Y',
                                                                      special_list='Y')
                star_phot['JD'] = jd
                target_list['JD'] = jd

                # dump the star phot file to a flux file
                star_phot[['TICID', 'ra', 'dec', 'x', 'y',
                           'JD', 'flux', 'flux_err', 'mag',
                           'mag_err', 'clean']].to_csv(Configuration.DIFFERENCED_DIRECTORY +
                                                       flx_file_nme, sep=" ", index=False)
                # dump the star phot file to a flux file
                target_list[['id', 'ra', 'dec', 'x', 'y', 'JD', 'flux', 'flux_err', 'clean']].\
                    to_csv(Configuration.DIFFERENCED_DIRECTORY +
                           file_prt[0] + '_' + Configuration.SPECIAL_LIST + '.flux', sep=" ", index=False)

        return

    @staticmethod
    def detrend_lightcurves(detrend_list):
        """ This function will detrend the stars in the raw light curve directory based on similarities in their
        shape, magnitude, and location.

        :parameter detrend_list - A data frame with the information of stars that are OK to detrend

        :return nothing is returned, but the detrended stars are output to the clean directory

        """

        # pull in the image quality data
        files = Utils.get_file_list(Configuration.CLEAN_DIRECTORY, '-' + Configuration.CAMERA +
                                    '-' + Configuration.CCD + '-' + Configuration.SECT_NUM + '-s_ffic_' +
                                    Configuration.FILE_EXT + 'd.fits')

        img_chk, ref_path = Master.check_frame_quality(files)
        img_chk['JD'] = img_chk.apply(lambda x: np.around(x['JD'], decimals=6), axis=1)

        # sort the detrend list based on magnitude
        detrend_stars = detrend_list.sort_values(by=['tessmag']).copy().reset_index(drop=True)

        # read in the light curves into a single data frame
        Utils.log("Detrending " + str(len(detrend_stars)) + " stars in sector: " + Configuration.SECTOR + " camera: " +
                  Configuration.CAMERA + " ccd: " + Configuration.CCD, "info", Configuration.LOG_SCREEN)

        # get the light curves for the initial number of light curves
        iter_num = 0
        raw_array = Photometry.read_in_lightcurves(Configuration.RAW_LC_DIRECTORY,
                                                   detrend_stars[0:Configuration.TREND_STARS]['TICID'].to_list())
        median_df = pd.DataFrame(raw_array - np.median(raw_array, axis=1).reshape(-1, 1))

        # make the gradient array
        grad_array = raw_array[:, 0:len(img_chk)-1] - raw_array[:, 1:len(img_chk)]
        grad_df = pd.DataFrame(grad_array)

        for idx, row, in detrend_stars.iterrows():

            # read in the raw light curve
            lc = pd.read_csv(Configuration.RAW_LC_DIRECTORY +
                             str(row.TICID.astype(int)) + "_" +
                             Configuration.SECTOR + "_" +
                             Configuration.CAMERA + "_" +
                             Configuration.CCD + ".lc",
                             names=['JD', 'clean', 'mag', 'err'], sep=" ")

            # get the light curve gradient
            grad_lc = lc['clean'][0:-1].to_numpy() - lc['clean'][1:].to_numpy()

            # get the spearman rank relationship for each light curve
            corr_lc = grad_df.apply(lambda x: x.corr(pd.Series(grad_lc), 'spearman'), axis=1)

            # check for like stars
            if len(corr_lc[(corr_lc > 0.85) & (corr_lc < .99999)]) > 50:
                # pull out stars with high spearman ranks correlation
                trend = np.average(median_df[(corr_lc > 0.85) & (corr_lc < 1)],
                                   weights=corr_lc[(corr_lc > 0.85) & (corr_lc < 1)], axis=0)
                lc['detrend'] = lc['clean'] - trend
            else:
                trend = np.average(median_df.iloc[np.argsort(corr_lc)[::-1][1:50]],
                                   weights=corr_lc.iloc[np.argsort(corr_lc)[::-1][1:50]], axis=0)
                lc['detrend'] = lc['clean'] - trend

                # remove any old versions of the light curve
                if os.path.exists(Configuration.DETREND_LC_DIRECTORY + str(detrend_stars.TICID.astype(int)) + "_" +
                                  str(Configuration.SECTOR) + "_" + str(Configuration.CAMERA) + "_" +
                                  str(Configuration.CCD) + ".lc"):
                    os.system('rm -f ' + Configuration.DETREND_LC_DIRECTORY + str(detrend_stars.TICID.astype(int)) +
                              "_" + str(Configuration.SECTOR) + "_" + str(Configuration.CAMERA) + "_" +
                              str(Configuration.CCD) + ".lc")

                # remove NaN values
                lc = lc.fillna(-9.999999)

                # write the cleaned light curve
                lc[['JD', 'detrend', 'clean', 'mag', 'err']].to_csv(Configuration.DETREND_LC_DIRECTORY +
                                                                    str(row.TICID.astype(int)) + "_" +
                                                                    str(Configuration.SECTOR) + "_" +
                                                                    str(Configuration.CAMERA) + "_" +
                                                                    str(Configuration.CCD) + ".lc", sep=" ",
                                                                    float_format='%.6f', header=False, index=False)

            # we are looking for stars within +/- 200 nearest magnitudes
            if (idx % 100 == 0) & (idx > 300):

                Utils.log("Moving the trend star filter by 100 stars.", "info", "Y")

                # set the bounds on the light curves to read in
                lw_bnd = Configuration.TREND_STARS + (iter_num * 100)
                up_bnd = Configuration.TREND_STARS + (iter_num * 100) + 100

                # get the light curves for the initial number of light curves
                raw_hold = Photometry.read_in_lightcurves(Configuration.RAW_LC_DIRECTORY,
                                                          detrend_stars[lw_bnd:up_bnd]['TICID'].to_list())
                median_hold = pd.DataFrame(raw_hold - np.median(raw_hold, axis=1).reshape(-1, 1))

                # add the next 100 stars, clip the top 100 and reorganize the data frame
                median_df = median_df.append(median_hold)
                median_df = median_df.drop(np.arange(100)).reset_index(drop=True)

                # make the gradient array
                grad_hold = raw_hold[:, 0:len(img_chk) - 1] - raw_hold[:, 1:len(img_chk)]
                grad_df = grad_df.append(pd.DataFrame(grad_hold))
                grad_df = grad_df.drop(np.arange(100)).reset_index(drop=True)
                iter_num += 1

            if (idx % 100 == 0) & (idx > 0):
                Utils.log("1000 stars have been de-trended in sector: " + Configuration.SECTOR + " camera: " +
                          Configuration.CAMERA + " ccd: " + Configuration.CCD + ". Now working on the next 1000. " +
                          str(len(detrend_stars) - idx - 1) + " stars remain to be cleaned.", "info",
                          Configuration.LOG_SCREEN)
        return

    @staticmethod
    def get_filtergraph_stars(master_list):
        """ This function will determine the list of stars we will want to detrend, primarily based on the expected
        TESSMAG, and the true tessmag.

        :parameter master_list - The star list from the master frame

        :return filtergraph_list - A list of stars to post to filtergraph
        """

        if os.path.exists(Configuration.MASTER_DIRECTORY + Configuration.SECTOR +
                          "_" + Configuration.CAMERA + "_" + Configuration.CCD + '_filtergraph_list.csv') == 1:
            Utils.log("Legacy filtergraph list found for Sector: " + Configuration.SECTOR +
                      " Camera: " + Configuration.CAMERA + " CCD: " + Configuration.CCD,
                      'info', Configuration.LOG_SCREEN)

            filtergraph_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.SECTOR +
                                           "_" + Configuration.CAMERA + "_" + Configuration.CCD +
                                           '_filtergraph_list.csv', index_col=0)

        else:
            Utils.log("Working to get the filtergraph list for Sector: " + Configuration.SECTOR +
                      " Camera: " + Configuration.CAMERA + " CCD: " + Configuration.CCD,
                      'info', Configuration.LOG_SCREEN)

            if Configuration.STAR_CHECK != -1:
                star_list = master_list[0:Configuration.STAR_CHECK].copy().reset_index(drop=True)
            else:
                star_list = master_list.copy().reset_index(drop=True)
                
            # split the data frame based on the zeropoint
            stars_zpt = star_list[star_list['tessmag'] <= Configuration.MG_CLP].copy().reset_index(drop=True)

            # calculate the expected zeropoint
            zpt = np.median(stars_zpt.tessmag.to_numpy() - stars_zpt.mag.to_numpy())
            zpt_std = np.std(stars_zpt.tessmag.to_numpy() - stars_zpt.mag.to_numpy())

            # determine which stars are OK to check for contamination
            stars_cratio = star_list.copy().reset_index(drop=True)
            stars_cratio['mag_offset'] = np.abs(stars_cratio.tessmag.to_numpy() -
                                                stars_cratio.mag.to_numpy() -
                                                zpt) / zpt_std

            # remove stars with wild zeropoint offsets
            filtergraph_list = stars_cratio[(stars_cratio['mag_offset'] < 2.5)].copy().reset_index(drop=True)
            filtergraph_list = filtergraph_list.sort_values('mag').reset_index(drop=True)

            # identify stars which slipped through the cracks during photometry
            index_to_drop = list()
            for idx, row in filtergraph_list.iterrows():
                # if the light curve doesn't exist, remove it
                if os.path.isfile(Configuration.RAW_LC_DIRECTORY + str(int(row['TICID'])) + "_" +
                                  Configuration.SECTOR + "_" + Configuration.CAMERA + "_" +
                                  Configuration.CCD + ".lc") is False:
                    index_to_drop.append(idx)

            # drop the non-existent light curves
            filtergraph_list = filtergraph_list.drop(index=index_to_drop).reset_index(drop=True)

            # write out the list of stars for filtergraph
            filtergraph_list[['TICID', 'x', 'y', 'mag', 'mag_err', 'tessmag']].to_csv(Configuration.MASTER_DIRECTORY +
                                                                                      Configuration.SECTOR + "_" +
                                                                                      Configuration.CAMERA + "_" +
                                                                                      Configuration.CCD +
                                                                                      '_filtergraph_list.csv')

            Utils.log(str(len(filtergraph_list)) + " stars found to post to filtergraph.",
                      'info', Configuration.LOG_SCREEN)

        return filtergraph_list

    @staticmethod
    def mk_raw_lightcurves(star_list):
        """ This function will create the individual raw light curve files for each star in the specific star list

        :parameter star_list - A data frame with the master frame flux data

        :return nothing is returned, but each light curve is output
        """
        if Configuration.SPECIAL_LIST != 'none':
            spec_files = Utils.get_file_list(Configuration.DIFFERENCED_DIRECTORY, '-' + Configuration.SECT + '-' +
                                        Configuration.CAMERA + '-' + Configuration.CCD + '-' +
                                        Configuration.SECT_NUM + '-s_ffic_' + Configuration.FILE_EXT + 'd' +
                                        '_' + Configuration.SPECIAL_LIST + '.flux')

            # combine the flux from the flux files, and write the raw light curves
            LightCurves.combine_special_list_files(Configuration.DIFFERENCED_DIRECTORY, spec_files)

        # get the file list of the differenced image flux information
        files = Utils.get_file_list(Configuration.DIFFERENCED_DIRECTORY, '-' + Configuration.SECT + '-' +
                                    Configuration.CAMERA + '-' + Configuration.CCD + '-' +
                                    Configuration.SECT_NUM + '-s_ffic_' + Configuration.FILE_EXT + 'd.flux')

        # combine the flux from the flux files, and write the raw light curves
        LightCurves.combine_flux_files(Configuration.DIFFERENCED_DIRECTORY, files, star_list)

        return

    @staticmethod
    def combine_flux_files(directory, files, star_list):
        """ This function combines all of the flux files in a given directory into a single data frame.

        :parameter directory - The directory where the files are located
        :parameter files - A list of the files
        :parameter star_list - The flux information for the master frame

        :returns data_df - A large data frame with all of the stellar flux information
        """

        # grab the master frame data
        master_df = star_list[['TICID', 'flux', 'flux_err']].copy().set_index('TICID')

        # there will be a length of 100,000 rows to read in at a time for a given flux file
        if Configuration.STAR_CHECK == -1:
            num_rrows = 100000
            nstars = len(star_list)
        else:
            num_rrows = 1000
            nstars = Configuration.STAR_CHECK

        # loop through the flux files, converting the flux to magnitude and then adding to the previous file
        for idx in range(0, nstars, num_rrows):

            Utils.log("Now working to create the raw light curves for Sector: " + Configuration.SECTOR +
                      " Camera: " + Configuration.CAMERA + " CCD:" + Configuration.CCD + ". Getting stars " +
                      str(idx) + " to " + str(num_rrows+idx) + ".", "info", Configuration.LOG_SCREEN)

            # update num_rrows if needed
            if nstars - idx < num_rrows:
                num_rrows = nstars - idx

            # make the holders for the light curves
            date = np.zeros(len(files))
            mags = np.zeros((len(files), num_rrows)) - 99.00
            clns = np.zeros((len(files), num_rrows)) - 99.00
            errs = np.zeros((len(files), num_rrows)) - 99.00

            for idy, file in enumerate(files):

                # read in the flux file
                if idx == 0:
                    flux_df = pd.read_csv(directory + file, sep=' ', index_col='TICID', header=0, nrows=num_rrows)
                if idx > 0:
                    flux_df = pd.read_csv(directory + file, sep=' ', index_col='TICID',
                                          header=0, skiprows=range(1, idx+1), nrows=num_rrows)

                # combine with the master frame photometry
                full_df = pd.merge(flux_df, master_df, on='TICID', how='left', suffixes=['', '_mast']).reset_index()

                # set the date
                date[idy] = full_df['JD'].iloc[0]

                # calculate the magnitudes
                mags[idy, :] = full_df['mag'].to_numpy()

                # calculate the clean magnitudes
                clns[idy, :] = full_df['clean'].to_numpy()

                # calculate the errors
                errs[idy, :] = full_df['mag_err'].to_numpy()

            # get the TICIDs for writing the light curves
            tics = full_df.TICID.to_numpy()

            # now write hte light curves
            LightCurves.write_light_curves(tics, date, mags, clns, errs)

        return


    @staticmethod
    def combine_special_list_files(directory, files):
        """ This function combines all of the flux files in a given directory into a single data frame.

        :parameter directory - The directory where the files are located
        :parameter files - A list of the files

        :returns data_df - A large data frame with all of the stellar flux information
        """

        target_list = pd.read_csv(Configuration.DATA_DIRECTORY + Configuration.SPECIAL_LIST + "_stars\\" +
                                  Configuration.SPECIAL_LIST + "_" + Configuration.SECTOR + ".csv", delimiter=',')

        target_list = target_list[(target_list.camera == int(Configuration.CAMERA)) &
                                  (target_list.ccd == int(Configuration.CCD))].copy().reset_index(drop=True)

        nstars = len(target_list)

        Utils.log("Now working to create the raw light curves for the Special List: " + Configuration.SPECIAL_LIST +
                  " in Sector: " + Configuration.SECTOR +
                  " Camera: " + Configuration.CAMERA + " CCD:" + Configuration.CCD + ".",
                  "info", Configuration.LOG_SCREEN)

        # make the holders for the light curves
        date = np.zeros(len(files))
        flxs = np.zeros((len(files), nstars))
        clns = np.zeros((len(files), nstars))
        errs = np.zeros((len(files), nstars))

        for idy, file in enumerate(files):

            # read in the flux file
            flux_df = pd.read_csv(directory + file, sep=' ', header=0)

            # set the date
            date[idy] = flux_df['JD'].iloc[0]

            # calculate the magnitudes
            flxs[idy, :] = flux_df['flux'].to_numpy()

            # calculate the clean magnitudes
            clns[idy, :] = flux_df['clean'].to_numpy()

            # calculate the errors
            errs[idy, :] = flux_df['flux_err'].to_numpy()

        # get the TICIDs for writing the light curves
        tics = flux_df['id'].to_numpy()

        # now write hte light curves
        LightCurves.write_light_curves(tics, date, flxs, clns, errs)

        return

    @staticmethod
    def write_light_curves(tics, date, mags, clns, errs):
        """ This function will write the times, magnitudes and errors to files.

        :parameter tics - The TICIDs for the current set of stars
        :parameter date - The dates for the light curves
        :parameter mags - The magnitudes for the light curves
        :parameter clns - The cleaned magnitude for the light curve
        :parameter errs - The errors in the light curves.

        :return - Nothing is returned, but the light curve files are written
        """

        # initialize the light curve data frame
        lc = pd.DataFrame(columns={'JD', 'clean', 'mag', 'err'})

        for idx, ticid in enumerate(tics):

            if (idx % 10000 == 0) and (idx > 0):
                Utils.log("10,000 light curves written for sector: " + str(Configuration.SECTOR) +
                          " Camera: " + str(Configuration.CAMERA) + " CCD: " + str(Configuration.CCD) +
                          ". Writing the next 10,000, " + str(len(tics) - idx - 1) +
                          " light curves remain to be written.", "info", Configuration.LOG_SCREEN)

            # add the time, magnitude and error to the data frame
            lc['JD'] = np.around(date, decimals=6)
            lc['clean'] = np.around(clns[:, idx], decimals=6)
            lc['mag'] = np.around(mags[:, idx], decimals=6)
            lc['err'] = np.around(errs[:, idx], decimals=6)

            # write the data to a text file
            if os.path.exists(Configuration.RAW_LC_DIRECTORY + str(ticid) + "_" + str(Configuration.SECTOR) + "_" +
                              str(Configuration.CAMERA) + "_" + str(Configuration.CCD) + ".lc"):
                os.system('rm -f ' + Configuration.RAW_LC_DIRECTORY + str(ticid) + "_" + str(Configuration.SECTOR) +
                          "_" + str(Configuration.CAMERA) + "_" + str(Configuration.CCD) + ".lc")

            # write the new file
            lc[['JD', 'clean', 'mag', 'err']].to_csv(Configuration.RAW_LC_DIRECTORY + str(ticid) + "_" +
                                                     str(Configuration.SECTOR) + "_" +
                                                     str(Configuration.CAMERA) + "_" +
                                                     str(Configuration.CCD) + ".lc",
                                                     sep=" ", header=False, index=False, na_rep='9.999999')

        return
