""" This class is used for basic functions not spcecific to this code base such as:
 logging, file writing, and testing."""
from config import Configuration
from libraries.utils import Utils
import psycopg2
import pandas as pd
from astroquery.mast import Catalogs
from astropy.coordinates import SkyCoord
import sys
import os
import json
import requests
from urllib.parse import quote as urlencode


class DBaccess:

    @staticmethod
    def get_ticid_query(path, ticids):
        """ This function will query the tic for all columns based on the TICID.

        :parameter path - The location of the SQL file
        :parameter ticids - The ticid of the stars

        :return A string with the SQL
        """

        # read in teh appropriate sql file
        sql = open(path, 'r')
        sql_cmd = sql.read()
        sql.close()

        # replace the necessary strings
        tics = ', '.join(ticids)
        sql_cmd = sql_cmd.replace('%(ticid)s', tics)

        return sql_cmd

    @staticmethod
    def get_starlist_query(path, cen_ra, cen_de, mx_dist):
        """ This function will replace the determined centers and distance in the SQL FILE and then return a string
        useful for querying the data base.

        :parameter path - The location of the SQL file
        :parameter cen_ra - The center RA coordinate
        :parameter cen_de - The center Dec coordinate
        :parameter mx_dist - The maximum edge to center distance

        :return A string with the SQL
        """

        # read in the appropriate sql file
        sql = open(path, 'r')
        sql_cmd = sql.read()
        sql.close()

        # replace the necessary strings
        sql_cmd = sql_cmd.replace('%(cen_ra)s', str(cen_ra))
        sql_cmd = sql_cmd.replace('%(cen_de)s', str(cen_de))
        sql_cmd = sql_cmd.replace('%(mx_dist)s', str(mx_dist))

        return sql_cmd

    @staticmethod
    def query_tic7_bulk(sql_cmd, stassunlab):
        """ This function will query tic v 7 on tessdev, and return a data frame based on the query used. This function
        will confirm a dumped .csv file does not exist prior to accessing the database.

        :parameter sql_cmd - The sql query to use
        :parameter stassunlab - The current computer being used for the program

        :return df - a data frame with the query results
        """
        # set up the connection object based on whether you are on tessdev
        if stassunlab != 'tessdev':
            conn = psycopg2.connect(host="10.2.188.37", # host="129.59.141.168",
                                    port=5432,
                                    database="tessdb",
                                    user="tessuser",
                                    password="4-users")
        if stassunlab == 'tessdev':
            conn = psycopg2.connect(host="10.2.188.37", # host="129.59.141.168",
                                    port=5432,
                                    database="tessdb",
                                    user="tessuser",
                                    password="4-users")
        # set up the cursor object
        cur = conn.cursor()

        Utils.log("Querying TICv7 for the next " + str(Configuration.BULK_QUERY) + " stars.",
                  "info", Configuration.LOG_SCREEN)

        # generate the data frame with the queried results
        df = pd.read_sql_query(sql_cmd, conn)

        Utils.log("Query complete!", "info", Configuration.LOG_SCREEN)

        # shut it down
        cur.close()
        conn.close()

        return df

    @staticmethod
    def query_tic7(sql_cmd, out_path, file_name, stassunlab):
        """ This function will query tic v 7 on tessdev, and return a data frame based on the query used. This function
        will confirm a dumped .csv file does not exist prior to accessing the database.

        :parameter sql_cmd - The sql query to use
        :parameter out_path - The output for the file
        :parameter file_name - The desired filename
        :parameter stassunlab - The current computer being used for the program

        :return df - a data frame with the query results
        """
        # set up the connection object based on whether you are on tessdev
        if stassunlab != 'tessdev':
            conn = psycopg2.connect(host="10.2.188.37",  #129.59.141.168",
                                    port=5432,
                                    database="tessdb",
                                    user="tessuser",
                                    password="4-users")
        if stassunlab == 'tessdev':
            conn = psycopg2.connect(host="10.2.188.37",
                                    port=5432,
                                    database="tessdb",
                                    user="tessuser",
                                    password="4-users")
        # set up the cursor object
        cur = conn.cursor()

        if os.path.isfile(out_path + file_name) == 1:
            Utils.log("Legacy file found, not querying TICv7.", "info", Configuration.LOG_SCREEN)
            # read in from a file
            df = pd.read_csv(out_path + file_name, index_col=0)
            Utils.log("CSV read complete.", "info", Configuration.LOG_SCREEN)

        if os.path.isfile(out_path + file_name) == 0:
            Utils.log("Querying TICv7...", "info", Configuration.LOG_SCREEN)

            # generate the data frame with the queried results
            df = pd.read_sql_query(sql_cmd, conn)

            Utils.log("Query complete, dumping to .csv file " + out_path + file_name, "info", Configuration.LOG_SCREEN)
            # dump the file to csv
            df.to_csv(out_path + file_name)

            Utils.log("Dump complete.", "info", Configuration.LOG_SCREEN)

        # convert the pk name to a TICID
        df = df.rename(columns={'pk': 'TICID'})

        # shut it down
        cur.close()
        conn.close()

        return df

    @staticmethod
    def mast_query(request):
        """Perform a MAST query. Lifted from: https://mast.stsci.edu/api/v0/pyex.html#MastCatalogsFilteredTicPy

            Parameters
            ----------
            request (dictionary): The MAST request json object

            Returns head,content where head is the response HTTP headers, and content is the returned data"""

        # Base API url
        request_url = 'https://mast.stsci.edu/api/v0/invoke'

        # Grab Python Version
        version = ".".join(map(str, sys.version_info[:3]))

        # Create Http Header Variables
        headers = {"Content-type": "application/x-www-form-urlencoded",
                   "Accept": "text/plain",
                   "User-agent": "python-requests/" + version}

        # Encoding the request as a json string
        req_string = json.dumps(request)
        req_string = urlencode(req_string)

        # Perform the HTTP request
        resp = requests.post(request_url, data="request=" + req_string, headers=headers)

        # Pull out the headers and response content
        head = resp.headers
        content = resp.content.decode('utf-8')

        return head, content

    @staticmethod
    def set_filters(parameters):
        """ Filtering function lifted from MAST website:
        https://mast.stsci.edu/api/v0/pyex.html#MastCatalogsFilteredTicPy """
        return [{"paramName": p, "values": v} for p, v in parameters.items()]

    @staticmethod
    def query_mast(cen_ra, cen_dec, region_size, out_path, file_name):
        """ This function will query the tic on mast instaed of the local tessdev databsae

        :parameter cen_ra - The ra center region
        :parameter cen_dec - The dec center region
        :parameter region_size - The region size of the query
        :parameter out_path - The directory to dump the data file
        :parameter file_name - The desired filename

        :return df - a data frame with the query results
        """

        if os.path.isfile(out_path + file_name) == 1:
            Utils.log("Legacy file found, not querying TIC on MAST", "info", Configuration.LOG_SCREEN)
            # read in from a file
            df = pd.read_csv(out_path + file_name, index_col=0)
            Utils.log("CSV read complete.", "info", Configuration.LOG_SCREEN)

        if os.path.isfile(out_path + file_name) == 0:
            Utils.log("Querying TIC on MAST...", "info", Configuration.LOG_SCREEN)

            # convert the coordinates to degrees
            ccd_region = SkyCoord(cen_ra, cen_dec, unit=('deg', 'deg'))

            # use the mast query functions to get the stars in 2MASS
            filts = DBaccess.set_filters({
                "typeSrc": ["tmgaia2"]
            })

            request = {"service": "Mast.Catalogs.Filtered.Tic.Position.Rows",
                       "format": "json",
                       "params": {
                           "columns": "ID,Tmag,ra,dec",
                           "filters": filts,
                           "ra": cen_ra,
                           "dec": cen_dec,
                           "radius": region_size
                       }}

            headers, out_string = DBaccess.mast_query(request)
            out_data = json.loads(out_string)
            df = pd.json_normalize(out_data['data'])
            df = df.rename(columns={'ID': 'TICID'})
            df = df.rename(columns={'Tmag': 'tessmag'})

            Utils.log("Query complete, dumping to .csv file " + out_path + file_name, "info", Configuration.LOG_SCREEN)

            # dump the file to csv
            df.to_csv(out_path + file_name)

            Utils.log("Dump complete.", "info", Configuration.LOG_SCREEN)

        return df
