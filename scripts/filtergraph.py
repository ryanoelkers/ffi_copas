""" This script work to generate the Filtergraph portal for release. It will also prepare the necessary files: tarballs,
plots, variable stats, csv files etc."""
from libraries.utils import Utils
from libraries.dbaccess import DBaccess
from config import Configuration
from scripts.lightcurves import LightCurves
import os
import pandas as pd
import matplotlib.pyplot as plt
import logging
import numpy as np
logging.getLogger('matplotlib.font_manager').disabled = True
from astroquery.mast import Catalogs


class Filtergraph:

    @staticmethod
    def make_filtergraph_table(star_list):
        """ This program makes the final release files and filtergraph table for upload.

        :parameter star_list - The list of stars to determine the filtergraph portal for

        :return nothing is returned, but the filtergraph table, figures, and variability statistics are finalized
        """

        # get the stars 'worth' de-trending
        filtergraph_list = LightCurves.get_filtergraph_stars(star_list)

        # get the varstats for the light curves
        Filtergraph.get_varstats(filtergraph_list)

        # make the figures for the light curves
        Filtergraph.plot_lightcurves(filtergraph_list)

        # make the filtergraph file from the relevant TIC
        Filtergraph.pull_filtergraph_portal(filtergraph_list)

        # combine all of the files to make one large file
        Filtergraph.make_filtergraph_file(filtergraph_list)

        return

    @staticmethod
    def make_filtergraph_file(filtergraph_list):
        """ This function will combine the varstats, filtergraph_list, and TIC file to make a file to
        upload to filtergraph. This file needs to be combined will all other files in order to upload to
        filtergraph appropriately.

        :parameter filtergraph_list - The list of stars which were detrended, and have filtergraph information.

        :return No parameter is returned, however the final file is combined appropriately.
        """

        # pull in the varstats file
        varstats = pd.read_csv(Configuration.VARSTATS_SECTOR_DIRECTORY + Configuration.SECTOR + "_" +
                               Configuration.CAMERA + "_" + Configuration.CCD + "_varstats.txt", delimiter=r'\s+')

        # append appropriate columns
        varstats['ticid'] = varstats.apply(lambda x: np.int(x['#Name'].split('_')[0].split('/')[-1]), axis=1)
        varstats['Lightcurve'] = 'Data_File'
        varstats['Figure'] = 'Plot_Lightcurve'
        varstats['Simbad'] = 'Simbad_Link'
        varstats['Aladin'] = 'Aladin_Link'
        varstats['Sector'] = Configuration.SECTOR
        varstats['Camera'] = Configuration.CAMERA
        varstats['CCD'] = Configuration.CCD

        # pull in the TIC file
        tic = pd.read_csv(Configuration.RELEASE_SECTOR_DIRECTORY + Configuration.SECTOR + "_" +
                          Configuration.CAMERA + "_" + Configuration.CCD + "_tic7_data.csv", delimiter=',', index_col=0)

        # combine with filtergraph_list
        if Configuration.WHERE_TO_QUERY =='local':
            full_tmp = pd.merge(filtergraph_list, tic, left_on='TICID', right_on='ticid')
            full = pd.merge(full_tmp, varstats, on='ticid')
            full = full.drop(columns=['#Name', 'pk', 'ticlinkpk', 'ticid', 'twomass_extkey',
                                      'sdss_extkey', 'kic', 'rpmjdwarf', 'tmpk'])

            # rename columns
            full = full.rename(columns={'TICID': 'TIC_ID', 'x': 'X_Pixel', 'y': 'Y_Pixel', 'mag': 'Inst_Mag',
                                        'mag_err': 'APER_Error', 'tessmag': "TESS_MAG"})
        else:
            full_tmp = pd.merge(filtergraph_list, tic, left_on='TICID', right_on='ID')
            full = pd.merge(full_tmp, varstats, left_on='TICID', right_on='ticid')
            full = full.drop(columns=['#Name', 'ID', 'ticid'])

            # rename columns
            full = full.rename(columns={'TICID': 'TIC_ID', 'x': 'X_Pixel', 'y': 'Y_Pixel', 'mag': 'Inst_Mag',
                                        'mag_err': 'APER_Error', 'tessmag': "TESS_MAG"})
        # re-order columns
        cols = list(full)
        cols.insert(1, cols.pop(cols.index('Lightcurve')))
        cols.insert(2, cols.pop(cols.index('Figure')))
        cols.insert(3, cols.pop(cols.index('Simbad')))
        cols.insert(4, cols.pop(cols.index('Aladin')))
        cols.insert(5, cols.pop(cols.index('Sector')))
        cols.insert(6, cols.pop(cols.index('Camera')))
        cols.insert(7, cols.pop(cols.index('CCD')))
        full = full.loc[:, cols]

        # dump the csv file to the release directory
        if (Configuration.CAMERA == '1') & (Configuration.CCD == '1'):
            full.to_csv(Configuration.RELEASE_SECTOR_DIRECTORY + Configuration.SECTOR + '_' + Configuration.CAMERA +
                        '_' + Configuration.CCD + '_filtergraph.csv', sep=',', index=False)
        # we will CAT the files later, so only include the header on the first file
        else:
            full.to_csv(Configuration.RELEASE_SECTOR_DIRECTORY + Configuration.SECTOR + '_' + Configuration.CAMERA +
                        '_' + Configuration.CCD + '_filtergraph.csv', sep=',', index=False, header=False)

        return

    @staticmethod
    def pull_filtergraph_portal(filtergraph_list):
        """ This function will query the necessary TIC for the release information.

        :parameter filtergraph_list - The list of stars which need tic information

        :return No value is returned, however the TIC file is dumped
        """

        if os.path.exists(Configuration.RELEASE_SECTOR_DIRECTORY + Configuration.SECTOR + "_" +
                          Configuration.CAMERA + "_" + Configuration.CCD + "_tic7_data.csv") is False:

            for idx in range(0, len(filtergraph_list), Configuration.BULK_QUERY):

                # create the appropriate index values
                idy = idx + Configuration.BULK_QUERY
                if idy > len(filtergraph_list):
                    idy = len(filtergraph_list)

                # pull the specific TICIDs to query
                tics = filtergraph_list['TICID'].iloc[idx:idy].astype(str).tolist()

                if Configuration.WHERE_TO_QUERY == 'local':
                    # create the correct SQL command
                    sql_cmd = DBaccess.get_ticid_query(Configuration.QUERIES_DIRECTORY + 'filtergraph.sql', tics)

                    # now query the TIC at this given distance and center coordinate position
                    df_query = DBaccess.query_tic7_bulk(sql_cmd, Configuration.MACHINE)
                else:

                    # now query MAST in bulks of 1000
                    df_query = Catalogs.query_criteria(catalog="Tic", ID=tics).to_pandas().dropna(subset=['TWOMASS'])
                    df_query = df_query.drop(columns=df_query.columns[list(df_query.columns).index("priority") + 1:])

                # either make a new data frame or append depending on the index
                if idx == 0:
                    df = df_query
                else:
                    df = df.append(df_query).reset_index(drop=True)

            # dump the file to csv
            df.to_csv(Configuration.RELEASE_SECTOR_DIRECTORY + Configuration.SECTOR + "_" +
                      Configuration.CAMERA + "_" + Configuration.CCD + "_tic7_data.csv")
        else:
            Utils.log("TIC file found for sector: " + Configuration.SECTOR + " camera: " + Configuration.CAMERA +
                      " ccd: " + Configuration.CCD + " skipping for now.", "info", Configuration.LOG_SCREEN)

        return

    @staticmethod
    def plot_lightcurves(filtergraph_list):
        """ This function will plot each light curve for release.

        :parameter filtergraph_list - The list of stars which will be posted to filtergraph

        :return Nothing is returned, however all of the stars will have their light curves plotted.
        """

        Utils.log("Now plotting light curves for sector: " + Configuration.SECTOR + " camera: " +
                  Configuration.CAMERA + " ccd: " + Configuration.CCD, "info", Configuration.LOG_SCREEN)

        for idx, tic in enumerate(filtergraph_list['TICID']):
            # simplify the filename for easy debugging
            nme = str(tic) + "_" + Configuration.SECTOR + "_" + Configuration.CAMERA + "_" + Configuration.CCD + ".png"

            # check to make sure the light curve isn't already plotted
            if os.path.isfile(Configuration.PLOTS_SECTOR_DIRECTORY + nme) is False:

                if os.path.isfile(Configuration.RAW_LC_DIRECTORY + str(tic) + "_" + Configuration.SECTOR + "_" +
                                  Configuration.CAMERA + "_" + Configuration.CCD + ".lc"):
                    # read in the specific light curve
                    lc = pd.read_csv(Configuration.RAW_LC_DIRECTORY + str(tic) + "_" + Configuration.SECTOR + "_" +
                                     Configuration.CAMERA + "_" + Configuration.CCD + ".lc", sep=" ",
                                     names=['JD', 'clean', 'mag', 'mag_err'])

                    # set the figure size
                    fig = plt.figure(figsize=(10, 10))

                    # make a scatter plot
                    plt.scatter(lc.JD, lc.clean, marker='.', c='k')

                    # set the labels
                    plt.xlabel('JD [days]', fontsize=15)
                    plt.ylabel('De-trended Instrumental Magnitude', fontsize=15)
                    plt.title('TICID: ' + str(tic))
                    plt.ylim([lc.mag.median() + lc.mag.std()*2.5, lc.mag.median() - lc.mag.std()*2.5])

                    # save the file
                    plt.savefig(Configuration.PLOTS_SECTOR_DIRECTORY + nme)
                    plt.close(fig)

                if (idx > 0) & (idx % 1000 == 0):
                    Utils.log("1000 light curves plotted. " + str(len(filtergraph_list) - idx - 1) +
                              " light curves remain.", "info", Configuration.LOG_SCREEN)

        return

    @staticmethod
    def get_varstats(filtergraph_list):
        """ This program will use VARTOOLS to generate the variability statistics for each light curve.

        :parameter filtergraph_list - The list of stars which will be posted to filtergraph

        :return No value is returned, but the varstats file is output
        """

        if os.path.exists(Configuration.VARSTATS_SECTOR_DIRECTORY + Configuration.SECTOR + "_" +
                          Configuration.CAMERA + "_" + Configuration.CCD + "_varstats.txt") is False:

            # make the file list
            nme = [Configuration.RAW_LC_DIRECTORY + str(s) + "_" + Configuration.SECTOR + "_" +
                   Configuration.CAMERA + "_" + Configuration.CCD + ".lc" for s in filtergraph_list['TICID'].to_list()]

            # dump to a text file for use with vartools
            pd.DataFrame(nme).to_csv(Configuration.VARSTATS_SECTOR_DIRECTORY + Configuration.SECTOR + "_" +
                                     Configuration.CAMERA + "_" + Configuration.CCD + "_starlist.txt",
                                     header=False, sep=" ", index=False)

            # get the initial light curve in the grouping to make the dates file
            lc = pd.read_csv(nme[0], sep=" ", names=['JD', 'clean', 'mag', 'mag_err'])

            # print out the dates to use later
            lc['JD'].to_csv(Configuration.VARSTATS_SECTOR_DIRECTORY + Configuration.SECTOR + "_" +
                            Configuration.CAMERA + "_" + Configuration.CCD + "_dates.csv", header=False, sep=" ")

            # run VARTOOLS on the files
            Utils.log("Calculating the variability statistics for sector: " + Configuration.SECTOR + " Camera: " +
                      Configuration.CAMERA + " CCD: " + Configuration.CCD, "info", Configuration.LOG_SCREEN)

            # send the vartools statement to the command line
            os.system("vartools -header -l " + Configuration.VARSTATS_SECTOR_DIRECTORY + Configuration.SECTOR + "_" +
                      Configuration.CAMERA + "_" + Configuration.CCD + "_starlist.txt " +
                      "-clip 2.5 1 -rms -Jstet 0.0208333 " +
                      Configuration.VARSTATS_SECTOR_DIRECTORY + Configuration.SECTOR + "_" + Configuration.CAMERA +
                      "_" + Configuration.CCD + "_dates.csv " +
                      "-rmsbin 1 60.0 -LS 0.1 27.0 0.1 1 0 -BLS q 0.01 0.1 0.1 27.0 " +
                      "2000 200 0 1 0 0 0 -readformat 0 1 2 4 > " +
                      Configuration.VARSTATS_SECTOR_DIRECTORY + Configuration.SECTOR + "_" + Configuration.CAMERA +
                      "_" + Configuration.CCD + "_varstats.txt")
#            os.system("vartools -header -l " + Configuration.VARSTATS_SECTOR_DIRECTORY + Configuration.SECTOR + "_" +
#                      Configuration.CAMERA + "_" + Configuration.CCD + "_starlist.txt " +
#                      "-clip 2.5 1 -rms -Jstet 0.0208333 " +
#                      Configuration.VARSTATS_SECTOR_DIRECTORY + Configuration.SECTOR + "_" + Configuration.CAMERA +
#                      "_" + Configuration.CCD + "_dates.csv " +
#                      "-rmsbin 1 60.0 -LS 0.1 27.0 0.1 1 0 -BLS q 0.01 0.1 0.1 27.0 " +
#                      "2000 200 0 1 0 0 0 -inputlcformat t:1,mag:2,err:4 > " +
#                      Configuration.VARSTATS_SECTOR_DIRECTORY + Configuration.SECTOR + "_" + Configuration.CAMERA +
#                      "_" + Configuration.CCD + "_varstats.txt")

            Utils.log("Variability calculation is complete.", "info", Configuration.LOG_SCREEN)

        else:
            Utils.log("Variability statistics have been calculated. Skipping.", "info", Configuration.LOG_SCREEN)

        return
