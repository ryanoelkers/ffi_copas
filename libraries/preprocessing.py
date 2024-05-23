from config import Configuration
import numpy as np
import astroalign as aa
from astropy.nddata.utils import Cutout2D
from astropy.wcs import WCS
from astropy.io import fits
from scipy.interpolate import griddata


class Preprocessing:

    @staticmethod
    def sky_subtract(img, header, sky_write='N'):
        """  The function breaks the image into bxs x bxs sections to save memeory, then cleans each section
        the sections are then recombined and smoothed to remove the transitions. The residual image is then subtracted
        from the image and the header is updated appropriately.

        :parameter img - The image to be cleaned
        :parameter header - The header object to be updated
        :parameter sky_write - Y/N if you want to write the residual background for de-bugging

        :return img_sub, header - The cleaned image and updated header file
        """

        # use the sampling space to make the appropriate size vectors
        lop = 2 * Configuration.PIX

        # size holder for later
        sze = int((Configuration.AXS / Configuration.PIX) * (Configuration.AXS / Configuration.PIX) +
                  (Configuration.AXS / Configuration.PIX) + (Configuration.AXS / Configuration.PIX) + 1)

        # calculate the sky statistics
        # sky_mean, sky_median, sky_sig = astropy.stats.sigma_clipped_stats(img, sigma=2.5)
        sky_median = np.nanmedian(img)
        sky_sig = np.nanstd(img)

        # create holder arrays for good and bad pixels
        x = np.zeros(shape=sze)
        y = np.zeros(shape=sze)
        v = np.zeros(shape=sze)
        s = np.zeros(shape=sze)
        nd = int(0)

        # begin the sampling of the "local" sky value
        for jj in range(0, Configuration.AXS + Configuration.PIX, Configuration.PIX):
            for kk in range(0, Configuration.AXS + Configuration.PIX, Configuration.PIX):
                il = np.amax([jj - lop, 0])
                ih = np.amin([jj + lop, Configuration.AXS - 1])
                jl = np.amax([kk - lop, 0])
                jh = np.amin([kk + lop, Configuration.AXS - 1])
                c = img[jl:jh, il:ih]

                # select the median value with clipping
                # lsky_mean, lsky, ssky = astropy.stats.sigma_clipped_stats(c, sigma=2.5)
                lsky_mean = np.nanmean(c)
                lsky = np.nanmedian(c)
                ssky = np.nanstd(c)
                x[nd] = np.amin([jj, Configuration.AXS - 1])  # determine the pixel to input
                y[nd] = np.amin([kk, Configuration.AXS - 1])  # determine the pixel to input
                v[nd] = lsky  # median sky
                s[nd] = ssky  # sigma sky
                nd = nd + 1

        # now we want to remove any possible values which have bad sky values
        rj = np.argwhere(v <= 0)  # stuff to remove
        kp = np.argwhere(v > 0)  # stuff to keep

        if len(rj) > 0:

            # keep only the good points
            xgood = x[kp]
            ygood = y[kp]
            vgood = v[kp]

            for jj in range(0, len(rj[0])):
                # select the bad point
                xbad = x[rj[jj]]
                ybad = y[rj[jj]]

                # use the distance formula to get the closest points
                rd = np.sqrt((xgood - xbad) ** 2. + (ygood - ybad) ** 2.)

                # sort the radii
                pp = sorted(range(len(rd)), key=lambda k: rd[k])

                # use the closest 10 points to get a median
                vnear = vgood[pp[0:9]]
                ave = np.median(vnear)

                # insert the good value into the array
                v[rj[jj]] = ave

        # now we want to remove any possible values which have bad sigmas
        rj = np.argwhere(s >= 2 * sky_sig)
        kp = np.argwhere(s < 2 * sky_sig)

        if len(rj) > 0:
            # keep only the good points
            xgood = np.array(x[kp])
            ygood = np.array(y[kp])
            vgood = np.array(v[kp])

            for jj in range(0, len(rj)):
                # select the bad point
                xbad = x[rj[jj]]
                ybad = y[rj[jj]]

                # use the distance formula to get the closest points
                rd = np.sqrt((xgood - xbad) ** 2. + (ygood - ybad) ** 2.)

                # sort the radii
                pp = sorted(range(len(rd)), key=lambda k: rd[k])

                # use the closest 10 points to get a median
                vnear = vgood[pp[0:9]]
                ave = np.median(vnear)
                if np.isfinite(ave) == 0:
                    ave = np.median(v[np.isfinite(v)])

                # insert the good value into the array
                v[rj[jj]] = ave

        # set up a meshgrid to interpolate to
        xi = np.linspace(0, Configuration.AXS - 1, Configuration.AXS)
        yi = np.linspace(0, Configuration.AXS - 1, Configuration.AXS)
        xx, yy = np.meshgrid(xi, yi)

        # remove any nan of inf values
        if len(v[~np.isfinite(v)]) > 0:
            v[~np.isfinite(v)] = np.median(v[np.isfinite(v)])

        # now we interpolate to the rest of the image with a cubic interpolation
        res = griddata((x, y), v, (xx, yy), method='cubic')

        # subtract the sky gradient and add back the median background
        img_sub = img - res
        fin_img = img_sub + np.median(sky_median)

        # update the header
        header['sky_medn'] = np.median(sky_median)
        header['sky_sig'] = np.median(sky_sig)
        header['sky_sub'] = 'yes'

        # if desired, write out the sky background to the working directory
        if sky_write == 'Y':
            fits.writeto('sky_background.fits', res, overwrite=True)

        return fin_img, header

    @staticmethod
    def align_img(img, header, ref_path):
        """ This function will align to a refernce position
        :parameter img - The image to flatten
        :parameter header - The image header file
        :parameter ref_path - The path to the reference iamge

        :return align_img, header - The updated image and header """
        # get the zeroth image for registration

        # read in the image
        ref, ref_header = fits.getdata(ref_path, header=True)

        ref_header['CRPIX1'] = Configuration.CRPIX
        ref_header['NAXIS1'] = Configuration.AXS
        ref_header['NAXIS2'] = Configuration.AXS

        # align the image
        # align_img = hcongrid(img, header, ref_header)
        align_img, footprint = aa.register(img, ref)

        # update the header
        header['CTYPE1'] = ref_header['CTYPE1']
        header['CTYPE2'] = ref_header['CTYPE2']
        header['CRVAL1'] = ref_header['CRVAL1']
        header['CRVAL2'] = ref_header['CRVAL2']
        header['CRPIX1'] = ref_header['CRPIX1']
        header['CRPIX2'] = ref_header['CRPIX2']
        header['CD1_1'] = ref_header['CD1_1']
        header['CD1_2'] = ref_header['CD1_2']
        header['CD2_1'] = ref_header['CD2_1']
        header['CD2_2'] = ref_header['CD2_2']
        header['aligned'] = 'Y'

        return align_img, header

    @staticmethod
    def bias_subtract(img, header):
        """ This function will subtract a bias frame
        :parameter img - The image to flatten
        :parameter header - The image header file

        :return bias_sub, header - The updated image and header """
        # read in the bias frame and subtract
        bias = fits.getdata(Configuration.CALIBRATION_DIRECTORY + 'bias.fits')

        # subtract the bias from the image
        bias_sub = img - bias

        # update the header
        header['bias_sub'] = 'Y'

        return bias_sub, header

    @staticmethod
    def flat_divide(img, header):
        """ This function will divide a flat field.
        :parameter img - The image to flatten
        :parameter header - The image header file

        :return flat_div, header - The updated image and header """
        # read in the flat frame
        flat = fits.getdata(Configuration.CALIBRATION_DIRECTORY + 'flat.fits')

        # subtract the bias from the image
        flat_div = img / flat

        # update the header
        header['flat_div'] = 'Y'

        return flat_div, header

    @staticmethod
    def mk_nme(file, difference_image='N', sky_subtract='N', bias_subtract='N', flat_divide='N', alignment='N'):
        """ This function will create the appropriate name for the file based on while steps are taken.
        :argument file - The string with the file name
        :argument difference_image - Y/N if image subtraction occured
        :argument sky_subtract - Y/N if sky subtraction was taken
        :argument bias_subtract - Y/N if a bias is subtracted
        :argument flat_divide - Y/N if a flat field is divided
        :argument alignment - Y/N if the image is aligned

        :return file_name - A string with the new file name
        """
        # make a new name for the file based on which actions are taken
        nme_hld = file.split('.')

        # if everything is N then the file name is the original filename
        file_name = file

        # update the file name with a 'd' if at the differencing step
        if difference_image == 'Y':
            file_name = nme_hld[0] + 'd.fits'

        # otherwise...
        if difference_image == 'N':
            # update the name to be appropriate for what was done to the file
            # nothing occurs
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'N'):
                file_name = nme_hld[0] + '.fits'
            # bias only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'N'):
                file_name = nme_hld[0] + '_b.fits'
            # flat
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'N'):
                file_name = nme_hld[0] + '_f.fits'
            # align
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'Y') and (sky_subtract == 'N'):
                file_name = nme_hld[0] + '_a.fits'
            # sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'Y'):
                file_name = nme_hld[0] + '_s.fits'

            # bias and flat only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'N'):
                file_name = nme_hld[0] + '_bf.fits'
            # bias and align only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'Y') and (sky_subtract == 'N'):
                file_name = nme_hld[0] + '_ba.fits'
            # bias and sky_subtract only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'Y'):
                file_name = nme_hld[0] + '_bs.fits'
            # bias and flat and align only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'N'):
                file_name = nme_hld[0] + '_bfa.fits'

            # flat and align only
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'N'):
                file_name = nme_hld[0] + '_fa.fits'
            # flat and sky_subtract only
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'Y'):
                file_name = nme_hld[0] + '_fs.fits'
            # flat and align and sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'Y'):
                file_name = nme_hld[0] + '_fsa.fits'

            # align and sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'Y') and (sky_subtract == 'Y'):
                file_name = nme_hld[0] + '_sa.fits'

            # all steps taken
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'Y'):
                file_name = nme_hld[0] + '_bfsa.fits'

        return file_name

    @staticmethod
    def remove_overscan(img, header):
        """ This function will remove the over scan region from the images
        :argument img - The np.array with the image
        :argument header - The header object

        :return cut_img, header - Return the new header and cut image"""

        # pull out the wcs in the header
        w = WCS(header)

        # cut the image, but maintain the appropriate header information and positioning
        cut = Cutout2D(img, (Configuration.XCUT, Configuration.YCUT), (Configuration.AXS, Configuration.AXS), wcs=w)

        # create a new image based on the cut data
        cut_img = cut.data

        # update the header
        header['CRPIX1'] = Configuration.CRPIX
        header['NAXIS1'] = Configuration.AXS
        header['NAXIS2'] = Configuration.AXS

        return cut_img, header
