# ===========================================================================================
# AUTHOR: M. Demuzere
# Date: 2021-07-15
# PURPOSE: Run TERRA-URB according to Urban Plumber requirements
# ===========================================================================================

# The script provides a generic way to run TERRA-URB for an urban site of interest.
# There are a number of different section, which should be checked carefully.
#	SECTION I:   simulation settings and folder locations
#	SECTION II:  Create the meteo input
#	SECTION III: Create param input files
# 	SECTION IV:  Run TERRA-URB
#	SECTION V:   Post-processing of outputfiles

# ===========================================================================================

SCRIPT_DIR = "/home/demuzmp4/Nextcloud/scripts/up-terra"
DATA_DIR   = "/home/demuzmp4/Nextcloud/data/up-terra"

import sys, os
import pandas as pd
import argparse
from argparse import RawTextHelpFormatter
# Import own modules. Reload done via auto reload (ipython required).
# Info: https://bit.ly/3xFGjBJ
sys.path.append(SCRIPT_DIR)
import utils as utils

parser = argparse.ArgumentParser(
    description="PURPOSE: Run TERRA for Urban Plumber Sites\n \n"
                "OUTPUT:\n"
                "- config: all model configuration details\n"
                "- input: TERRA readable meteo forcing\n"
                "- output: all results",
    formatter_class=RawTextHelpFormatter
)

# Required arguments
parser.add_argument(type=str, dest='SITE',
                    help='Name of Urban Plumber site',
                    )
parser.add_argument(type=str, dest='EXP',
                    help="Name of Experiment: \n"
                         "* 'b': BASELINE\n"
                         "* 'd': DETAILED\n"
                         "* 'lcz': optional, to run with LCZ-based parameters\n",
                    )
parser.add_argument(type=str, dest='VERSION',
                    help='Identifier to keep track of experiment version',
                    )
parser.add_argument(type=str, dest='START_RUN',
                    help="Defines where to start the run. Options: \n"
                         "* 'coverage': start at beginning, including all spin-up\n"
                         "* 'analysis': run analysis period only\n"
                         "* 'analysis_Xyr': run X years of spin-up before analysis\n"
                         "* 'analysis_Xyr_production': as analysis_Xyr, but filling all timestamps to time_coverage_start\n"
                         "* 'debug': run only last 90 days of analysis period\n",
                    )

args = parser.parse_args()

# Arguments to script
SITE = args.SITE
EXP = args.EXP
VERSION = args.VERSION
START_RUN = args.START_RUN

# Arguments for testing
# SITE = 'SG-TelokKurau'
# SITE = 'AU-Preston'
# EXP = 'd'
# VERSION = 'wpc'
# START_RUN = 'debug'
# CASE = 'urb'

# Fixed arguments
TERRA_VERSION  =  "TERRA_4.11"
MET_FORC_DIR_V = "Urban-PLUMBER_sites_nc_v1-1"

for CASE in ['urb','veg']:

    print("=========================================")
    print(f"START simulation: {SITE}, {CASE}, {EXP}, {VERSION}")
    print("=========================================")


    # ************************************************************************
    #		SECTION I: SIMULATIONS AND FOLDERS
    # ************************************************************************

    # Read site look up: Contains LCZ, MET_FORC_V
    site_lookup = pd.read_csv(f"{SCRIPT_DIR}/tables/site_lookup.csv", index_col=0)

    # Create info dir that contains all relevant paths and data.
    run_info = utils.set_run_info(
        SITE=SITE,
        CASE=CASE,
        EXP=EXP,
        VERSION=VERSION,
        SCRIPT_DIR=SCRIPT_DIR,
        DATA_DIR=DATA_DIR,
        TERRA_VERSION=TERRA_VERSION,
        LCZ=site_lookup.loc[SITE,'LCZ'],
        MET_FORC_DIR_V=MET_FORC_DIR_V,
        MET_FORC_V=site_lookup.loc[SITE,'MET_FORC_V'],
        START_RUN= START_RUN,
        TO_LOCAL_TIME=False,
        OUTPUT_TYPE='GRB', # 'GRB' or 'ASC'
        DISP_ROUGH='' ,# '' or '_mac' or '_kanda'
        )

    for i_dir in [run_info['SIM_CFG'],
                  run_info['OUTPUTDIR'],
                  run_info['SIM_INPUT']]:
        if not os.path.exists(i_dir):
            os.makedirs(i_dir)

    if not os.path.exists(run_info['RUN_INFO']):
        with open(run_info['RUN_INFO'], 'w') as f:
            for key in run_info.keys():
                f.write("%s,%s\n"%(key,run_info[key]))

    #************************************************************************
    # SECTION II: CREATE METEO INPUT (.bin FILE)
    #************************************************************************

    # Always create input for experiment, to make sure it is correct
    utils.create_forcing_bin(
        run_info=run_info)


    #************************************************************************
    # SECTION III: CREATE PARAM INPUT FILES
    #************************************************************************
    utils.create_params(
        run_info=run_info
    )

    utils.put_input_params(
        run_info=run_info
    )

    #************************************************************************
    # SECTION IV: RUN THE MODEL
    #************************************************************************
    utils.compile_model(
        run_info=run_info
    )

    utils.run_model(
        run_info=run_info
    )

    #************************************************************************
    # SECTION V: POST-PROCESSING OUTPUT
    #************************************************************************
    output_file = utils.convert_to_netcdf(
        run_info=run_info)

    print("===========================================================")
    print(f"Done for {CASE}. \n"
          f"Output file: {output_file}")
    print("===========================================================")

print("\n")
print("===========================================================")
print(f"Post-processing started: {SITE}, {EXP}, {VERSION}")
print("===========================================================")

# # Fixed parameters
soil_depths = [0.01, 0.03, 0.09, 0.27, 0.81, 2.43, 7.29, 21.87]
num_soil_layers = len(soil_depths)
missing_float   = -9999.

# Read run_info file.
run_info = pd.read_csv(os.path.join(
    f"{DATA_DIR}/output/{SITE}/output",
    VERSION,
    EXP,
    'run_info.csv'
),index_col=0, header=None, squeeze=True).to_dict()

# Read the meta_data file.
sim_info = pd.read_csv(run_info['SIM_INFO'],index_col=0, squeeze=True).to_dict()

# Get dataframe with results + forcing (GRBOUT)
data, data_SoilMoist, data_SoilTemp = utils.create_final_df_from_grb(
    run_info=run_info,
    sim_info=sim_info)

# Create empty netcdf to store results in
utils.create_empty_netcdf(
    run_info= run_info,
    sim_info= sim_info,
    data=data,
    num_soil_layers=num_soil_layers,
    missing_float=missing_float)

# Put results in this netcdf
utils.set_netcdf_data(
    run_info=run_info,
    data=data,
    data_SoilMoist=data_SoilMoist,
    data_SoilTemp=data_SoilTemp,
    num_soil_layers=num_soil_layers,
    missing_float=missing_float)

print("===========================================================")
print(f"Post-processing done, following Urban-Plumber requirements\n"
      f"Final file: {run_info['FINAL_FNAME_OUT']}")
print("===========================================================")
