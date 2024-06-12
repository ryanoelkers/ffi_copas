""" This serves as the configuration file for the FFI reduction. The goal is to only change these values when it is
necessary to update for a reduction."""


class Configuration:
    # Debugging loops
    STAR_CHECK = 25000  # this is the number of stars to run photometry on the frame SET TO -1 UNLESS DEBUGGING

    # Computer for reduction
    MACHINE = 'PHY-MASS206'
    SYSTEM = 'linux'

    BULK_QUERY = 25000  # 1000 for tessdb probably closer to 25,000 for MAST
    WHERE_TO_QUERY = 'MAST'
    SPECIAL_LIST = 'Y'
    SPECIAL_LIST_NAME = 'copas'

    # Sector, Camera, CCD information
    SECTOR = 'sector70'
    SECTOR_NUMBER = '0070'
    CAMERA = '4'
    CCD = '4'
    SECT = 's0070'  # 's0011' # 's0010'
    SECT_NUM = '0265'  # '0143'  # '0140'

    # logging information, i.e. print to the screen?
    LOG_SCREEN = 'Y'

    # do you want to sky subtract, bias subtract, flatten, or align the images?
    BIAS_SUBTRACT = "N"
    FLAT_DIVIDE = 'N'
    SKY_SUBTRACT = 'Y'
    ALIGNMENT = 'N'

    # which reduction step would you like to start at?
    CLEAN_SKIP = 'Y'
    MASTER_SKIP = 'N'
    DIFFERENCE_SKIP = 'N'
    PHOTOMETRY_SKIP = 'Y'
    MAKE_RAW_LIGHTCURVE_SKIP = 'Y'
    MAKE_DETREND_LIGHTCURVE_SKIP = 'Y'
    MAKE_FILTERGRAPH_PORTAL_SKIP = 'Y'
    MAST_RELEASE_SKIP = 'Y'

    # what type of master frame would you like (single or multi)
    MASTER_TYPE = 'single'

    # update image specific information (i.e. axis size - for TESS ~2048x2048)
    AXS = 2048
    CRPIX = 1001.
    XCUT = 1068
    YCUT = 1024

    # update sky subtraction specific information
    PIX = 32

    # update photometry specific information
    APER_SIZE = 3.5
    ANNULI_INNER = 5.
    ANNULI_OUTER = 7.
    ERROR_LIMIT = 0.5

    # update here for detrend star lists
    MX_SEP = APER_SIZE * 2
    MG_CLP = 13
    DMG_CLP = 3.
    CONTAM_LIMIT = 5.
    TREND_LIMIT = 10.
    TREND_STARS = 500
    TREND_SCALE_LIMIT = 1.25
    CORRELATION_SCORE = 0.85

    # update the differencing information, primarily the number of stars to use, and the kernel size
    KRNL = 2  # kernel size 2 * KNRL + 1
    STMP = 3  # stamp size ot use 2 * STMP + 1
    ORDR = 0  # order of the kernel to use, 0 is stationary, 1 or 2 is spatially varying
    NRSTARS = 1000  # number of stars used to solve for kernel
    BRIGHT_STARS = 20000  # the top stars to search for in kernel stars
    KERNEL_LIMIT = 0.5  # the maximum allowable offset in zeropoint in magnitudes
    AXS_LIMIT = 30  # the number of pixel close to the edge of the frame to use
    RMS_LOW_LIMIT = 0.005  # the lower limit on precision to use for the kernel stars
    RMS_UP_LIMIT = 0.02  # the upper limit on precision to use for the kernel stars

    # output paths for logging, temporary files, figures etc
    WORKING_DIRECTORY = '/home/oelkerrj/PycharmProjects/ffi_copas/'
    DATA_DIRECTORY = '/media/oelkerrj/DATA/copas/'
    ANALYSIS_DIRECTORY = WORKING_DIRECTORY + 'analysis/'
    LEGACY_DIRECTORY = WORKING_DIRECTORY + 'legacy/'
    LIBRARY_DIRECTORY = WORKING_DIRECTORY + 'libraries/'
    LOG_DIRECTORY = WORKING_DIRECTORY + 'logs/'
    CALIBRATION_DIRECTORY = WORKING_DIRECTORY + 'calibration/'
    QUERIES_DIRECTORY = WORKING_DIRECTORY + 'queries/'

    # output directories for data products
    CLEAN_DIRECTORY = DATA_DIRECTORY + SECTOR + '/clean/'
    RAW_DIRECTORY = DATA_DIRECTORY + SECTOR + '/raw/'
    DIFFERENCED_DIRECTORY = DATA_DIRECTORY + SECTOR + '/diff/'
    MASTER_DIRECTORY = DATA_DIRECTORY + SECTOR + '/master/'
    LC_DIRECTORY = DATA_DIRECTORY + SECTOR + '/lc/'
    RAW_LC_DIRECTORY = LC_DIRECTORY + 'raw/'
    RAW_LC_SPEC_DIRECTORY = RAW_LC_DIRECTORY + SPECIAL_LIST + '/'
    DETREND_LC_DIRECTORY = LC_DIRECTORY + 'detrend/'
    RELEASE_DIRECTORY = DATA_DIRECTORY + 'release/'
    RELEASE_SECTOR_DIRECTORY = RELEASE_DIRECTORY + SECTOR + '/'
    MAST_DIRECTORY = RELEASE_SECTOR_DIRECTORY + 'mast/'
    MAST_SECTOR_DIRECTORY = RELEASE_SECTOR_DIRECTORY + 's' + SECTOR_NUMBER + '/'
    PLOTS_DIRECTORY = DATA_DIRECTORY + 'plots/'
    PLOTS_SECTOR_DIRECTORY = PLOTS_DIRECTORY + SECTOR + '/'
    VARSTATS_DIRECTORY = DATA_DIRECTORY + 'varstats/'
    VARSTATS_SECTOR_DIRECTORY = VARSTATS_DIRECTORY + SECTOR + '/'

    # directories for coding
    CODE_DIFFERENCE_INIT_DIRECTORY = DATA_DIRECTORY + 'difference/'
    CODE_DIFFERENCE_INIT_SECT_DIRECTORY = DATA_DIRECTORY + 'difference/' + SECTOR + '/'
    CODE_DIFFERENCE_CAM_DIRECTORY = DATA_DIRECTORY + 'difference/' + SECTOR + '/cam' + CAMERA + '/'
    CODE_DIFFERENCE_DIRECTORY = DATA_DIRECTORY + 'difference/' + SECTOR + '/cam' + CAMERA + '/ccd' + CCD + '/'

    # directory_list
    DIRECTORIES = [ANALYSIS_DIRECTORY, DATA_DIRECTORY, LEGACY_DIRECTORY, LIBRARY_DIRECTORY, LOG_DIRECTORY,
                   CLEAN_DIRECTORY, RAW_DIRECTORY, DIFFERENCED_DIRECTORY, LC_DIRECTORY, RAW_LC_DIRECTORY,
                   DETREND_LC_DIRECTORY, MASTER_DIRECTORY, CALIBRATION_DIRECTORY, CODE_DIFFERENCE_INIT_DIRECTORY,
                   CODE_DIFFERENCE_INIT_SECT_DIRECTORY, CODE_DIFFERENCE_CAM_DIRECTORY, CODE_DIFFERENCE_DIRECTORY,
                   RELEASE_DIRECTORY, RELEASE_SECTOR_DIRECTORY, PLOTS_DIRECTORY, PLOTS_SECTOR_DIRECTORY,
                   VARSTATS_DIRECTORY, VARSTATS_SECTOR_DIRECTORY, QUERIES_DIRECTORY, RAW_LC_SPEC_DIRECTORY]

    # file extension
    FILE_EXT = ''
    if BIAS_SUBTRACT == 'Y':
        FILE_EXT = FILE_EXT + 'b'
    if FLAT_DIVIDE == 'Y':
        FILE_EXT = FILE_EXT + 'f'
    if SKY_SUBTRACT == 'Y':
        FILE_EXT = FILE_EXT + 's'
    if ALIGNMENT == 'Y':
        FILE_EXT = FILE_EXT + 'a'
