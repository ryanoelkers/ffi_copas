from scripts.clean import Clean
from scripts.master import Master
from scripts.difference import BigDiff
from scripts.lightcurves import LightCurves
from scripts.filtergraph import Filtergraph
from scripts.mast_release import MastRelease
from scripts.detection import Detection
from libraries.utils import Utils
from config import Configuration
import pandas as pd
import matplotlib
matplotlib.use('Agg')

# do necessary prep work
Utils.create_directories(Configuration.DIRECTORIES)

if Configuration.CLEAN_SKIP == 'N':
    # begin cleaning the images
    Clean.clean_images(sky_subtract=Configuration.SKY_SUBTRACT,
                       bias_subtract=Configuration.BIAS_SUBTRACT,
                       flat_divide=Configuration.FLAT_DIVIDE,
                       alignment=Configuration.ALIGNMENT)
else:
    Utils.log("Skipping image cleaning.", "info", Configuration.LOG_SCREEN)

if Configuration.MASTER_SKIP == 'N':
    # make the master frame and get the basic information
    star_list = Master.mk_master(Configuration.MASTER_TYPE)
else:
    # if you skip the master frame pull in the star list just in case
    star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.SECTOR + "_" +
                            Configuration.CAMERA + '_' + Configuration.CCD + '_master.data', index_col=0)
    Utils.log("Skipping master file data pull.", "info", Configuration.LOG_SCREEN)

if Configuration.DIFFERENCE_SKIP == 'N':
    # difference the images
    BigDiff.difference_images(star_list)
else:
    Utils.log("Skipping image differencing.", "info", Configuration.LOG_SCREEN)

if Configuration.PHOTOMETRY_SKIP == 'N':
    # get the flux photometry for the images
    LightCurves.diff_img_phot(star_list)
else:
    Utils.log("Skipping photometry of differenced images.", "info", Configuration.LOG_SCREEN)

if Configuration.MAKE_DETECTION_FRAME_SKIP == 'N':
    # make the detection frames for the given sector, camera, ccd
    Detection.make_detection_frames(Configuration.BIN_SCALE)
else:
    Utils.log("Skipping making detection frames.", "info", Configuration.LOG_SCREEN)

if Configuration.MAKE_RAW_LIGHTCURVE_SKIP == 'N':
    # create the raw light curves from the flux files
    LightCurves.mk_raw_lightcurves(star_list)
else:
    Utils.log("Skipping making raw light curves.", "info", Configuration.LOG_SCREEN)

if Configuration.MAKE_DETREND_LIGHTCURVE_SKIP == 'N':

    # get the variable statistics for all stars
    LightCurves.detrend_lightcurves(star_list)
else:
    Utils.log("Skipping making de-trended light curves.", "info", Configuration.LOG_SCREEN)

# generate the filtergraph tic / varstats file
if Configuration.MAKE_FILTERGRAPH_PORTAL_SKIP == 'N':

    # make the filtergraph table
    Filtergraph.make_filtergraph_table(star_list)

else:
    Utils.log("Skipping making the filtergraph portal.", "info", Configuration.LOG_SCREEN)

# generate the files for MAST release
if Configuration.MAST_RELEASE_SKIP == 'N':

    # generate the .fits files for release
    MastRelease.make_mast_lc_bulk(star_list)
else:
    Utils.log("Skipping making MAST release files.", "info", Configuration.LOG_SCREEN)
Utils.log("All done! See ya later, alligator.", "info", Configuration.LOG_SCREEN)
