import shutil
import numpy as np
import pandas as pd
import netCDF4 as nc
import h5py
import xarray as xr
import os, sys
import re
import time
import glob
from scipy.special import erf
from datetime import timedelta

def set_run_info(SITE: str,
                 CASE: str,
                 EXP: str,
                 VERSION: str,
                 SCRIPT_DIR: str,
                 DATA_DIR: str,
                 TERRA_VERSION: str,
                 LCZ: int,
                 MET_FORC_DIR_V: str,
                 MET_FORC_V: str,
                 START_RUN: str,
                 TO_LOCAL_TIME=False,
                 OUTPUT_TYPE='GRB',
                 DISP_ROUGH='',
                 ) -> dict:

    # Initialize the run_info dictionary
    run_info = {}

    # Set INPUT/OUTPUT interval as retrieved from netcdf file
    FORCING_FILE = f"{DATA_DIR}/input/{MET_FORC_DIR_V}/{SITE}_metforcing_{MET_FORC_V}.nc"
    ds = xr.open_dataset(FORCING_FILE)
    run_info['INPUT_INTERVAL']  = int(ds.attrs['timestep_interval_seconds']/60)
    run_info['OUTPUT_INTERVAL'] = run_info['INPUT_INTERVAL']

    run_info.update({
        'TERRA_CFG':      f"{SCRIPT_DIR}/terra/src_config",
        'SRCDIR':         f"{SCRIPT_DIR}/terra/{TERRA_VERSION}/src",
        'INPUT_NAMELIST': f"{SCRIPT_DIR}/terra/src_config/input.namelist.temp",
        'INPUT_ARRAY':    f"{SCRIPT_DIR}/terra/src_config/INPUT_TERRA.temp",
        'INPUT_DATA_SOIL':f"{SCRIPT_DIR}/terra/src_config/data_soil-temp.f90",

        'MET_FORCING':    f"{DATA_DIR}/input/{MET_FORC_DIR_V}/{SITE}_metforcing_{MET_FORC_V}.nc",

        'SIM_CFG':        f"{DATA_DIR}/output/{SITE}/config/{VERSION}/{EXP}/{CASE}",
        'SIM_INPUT':      f"{DATA_DIR}/output/{SITE}/input/{VERSION}/{EXP}/{run_info['INPUT_INTERVAL']}",
        'OUTPUTDIR':      f"{DATA_DIR}/output/{SITE}/output/{VERSION}/{EXP}/{CASE}",
        'SIM_INFO':       f"{DATA_DIR}/output/{SITE}/output/{VERSION}/{EXP}/sim_info.csv",
        'RUN_INFO':       f"{DATA_DIR}/output/{SITE}/output/{VERSION}/{EXP}/run_info.csv",
        'YUDRDAT':        f"{DATA_DIR}/output/{SITE}/config/{VERSION}/{EXP}/{CASE}/YUDRDAT",

        'LOOKUPDIR':      f"{SCRIPT_DIR}/tables",
        'UCP_FILE':       f"{SCRIPT_DIR}/tables/LCZ_UCP_default.csv",
        'AHF_FILE':       f"{SCRIPT_DIR}/tables/AHF_2005_2.5min.nc",

        'TO_LOCAL_TIME':  TO_LOCAL_TIME,
        'START_RUN':      START_RUN,
        'SITE':           SITE,
        'VERSION':        VERSION,
        'EXP':            EXP,
        'CASE':           CASE,
        'OUTPUT_TYPE':    OUTPUT_TYPE,
        'LCZ':            int(LCZ),
        'DISP_ROUGH':     DISP_ROUGH,
        'TERRA_VERSION':  TERRA_VERSION,
        'FINAL_FNAME_OUT':f"{DATA_DIR}/output/{SITE}/output/{TERRA_VERSION}_{SITE}_{EXP}_{VERSION}.nc"

    })

    if len(run_info['OUTPUTDIR']) > 88:
        print("WARNING: output directory name too long.\n"
              "Execution of terra_exec will fail.\n\n"
              "Try using a shorter VERSION name, eg. DB1 (in case of debug)\n\n"
              "Exiting now ...")
        sys.exit()

    # Define starting date of run here.
    ds = xr.open_dataset(run_info['MET_FORCING'])

    # Full run
    if START_RUN == 'coverage':
        run_info.update({'START_DATE':
                             pd.to_datetime(ds.attrs['time_coverage_start'])})
    # Only analysis period, no spin-up
    elif START_RUN == 'analysis':
        startdate = ds.attrs['time_analysis_start']
        startdate = pd.to_datetime(startdate).replace(hour=0, minute=0)
        run_info.update({'START_DATE': startdate})
    # Analysis period with X year(s) of spin up
    elif 'yr' in START_RUN:
        if not '_production' in START_RUN:
            spinup_years = int(START_RUN.replace('analysis_','').replace('yr',''))
        else:
            spinup_years = int(START_RUN.replace('analysis_', '')
                               .replace('yr', '')
                               .replace('_production', ''))
        startdate = ds.attrs['time_analysis_start']
        startdate = pd.to_datetime(startdate).replace(hour=0, minute=0)
        startdate = startdate - timedelta(days=spinup_years*365)
        run_info.update({'START_DATE': startdate})
    # Last 90 days of analysis, for debugging.
    elif START_RUN == 'debug':
        startdate = ds.attrs['time_coverage_end']
        startdate = pd.to_datetime(startdate) - timedelta(days=90)
        startdate = startdate.replace(hour=0, minute=0)
        run_info.update({'START_DATE': startdate})

    print(f">>>> RUN STARTING AT: {run_info['START_DATE']}\n")
    return run_info


def _create_forcing(run_info):
    '''
    Parameters
    ----------
    info (dictionary): script information

    Output
    ------
    forcing (dataframe)  : forcing data in model input form
    formats (dictionary) : string formatters for each dataframe column
    '''

    # load forcing information from netcdf
    ds = xr.open_dataset(run_info['MET_FORCING'])

    # Start run based on predefined settings
    time_coverage_end         = ds.attrs['time_coverage_end']
    timestep_interval_seconds = ds.attrs['timestep_interval_seconds']
    times = pd.date_range( start = run_info['START_DATE'], # time_coverage_start
                           end   = time_coverage_end,
                           freq  = '%sS' %timestep_interval_seconds)

    # Reduce input data to selected times
    ds = ds.sel(time=slice(times[0], times[-1]))

    Wind_E = ds['Wind_E'].values.flatten()
    Wind_N = ds['Wind_N'].values.flatten()
    Tair   = ds['Tair'].values.flatten()
    Qair   = ds['Qair'].values.flatten()
    PSurf  = ds['PSurf'].values.flatten()
    SWdown = ds['SWdown'].values.flatten()
    LWdown = ds['LWdown'].values.flatten()
    Rainf  = ds['Rainf'].values.flatten()
    Snowf  = ds['Snowf'].values.flatten()

    # change to local time if model requires this
    # Not needed for TERRA_URB
    if eval(str(run_info['TO_LOCAL_TIME'])):
        print(">>>> Use of Local Time Enabled <<<<")
        local_utc_offset_hours = ds.attrs['local_utc_offset_hours']
        offset = pd.Timedelta('%s hours' %local_utc_offset_hours)
        times = times + offset

    # create empty dataframe
    forcing = pd.DataFrame(index=times)

    # derived
    wind = np.sqrt(Wind_E**2 + Wind_N**2)
    #prec = (Rainf+Snowf)*timestep_interval_seconds
    prec = (Rainf + Snowf)
    qd2m = Qair

    # fill dataframe with forcing data with appropriate conversions
    forcing = forcing.assign(
                        glob = SWdown,          # [W/m2]
                        trat = LWdown,          # [W/m2]
                        prec = prec*3600,       # [kg/m2/h]
                        prss = PSurf/100,       # [hPa]
                        qd2m = qd2m*1000,       # [g/kg], Otherwise rounding errors?
                        te2m = Tair,            # [K]
                        w10m = wind,            # [m/s]
                        )

    # set column formats
    # include comma (,) in string if csv required e.g. '{:04d},'
    # for more formatting info see https://pyformat.info
    formats = {
            'glob'  : '{:8.2f}'.format,
            'trat'  : '{:8.2f}'.format,
            'prec'  : '{:8.2f}'.format,
            'prss'  : '{:8.2f}'.format,
            'qd2m'  : '{:8.2f}'.format,
            'te2m'  : '{:8.2f}'.format,
            'w10m'  : '{:8.2f}'.format,
            }

    return forcing, formats


def create_forcing_bin(run_info):

    print('creating model forcing data')
    data_df, formats = _create_forcing(run_info)

    print('writing forcing files')
    for var in ['glob','prec','prss','qd2m','te2m','trat','w10m']:

        # If files already exist, remove
        dat_file = f"{run_info['SIM_INPUT']}/{var}.dat"
        if os.path.exists(dat_file):
            os.remove(dat_file)

        # create single-variable dataframe where each row = one day
        # df = pd.DataFrame([group[1].values for group in
        #                    data_df[var].groupby(data_df.index.dayofyear)])
        df = pd.DataFrame([group[1].values for group in
                           data_df[var].groupby(data_df.index.date)])
        # reset index to start at 1
        df.index = df.index + 1

        # Set index values as first columns
        df.insert(loc=0, column='index', value=df.index)

        # Get the number of days, see below
        nr_days = df.shape[0]

        # set nans to zero, at end of file. Should be clipped from analysis
        df = df.fillna(0)

        with open(dat_file, 'w') as file:
            file.writelines('\n') #Blank line needed
            #file.writelines(data_str)

            if int(run_info['INPUT_INTERVAL']) == 60:
                fmt = '%4d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f'\
                      '%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f'\
                      '%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f'
            elif int(run_info['INPUT_INTERVAL']) == 30:
                fmt = '%4d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f'\
                      '%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f'\
                      '%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f'\
                      '%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f'\
                      '%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f'\
                      '%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f'
            np.savetxt(file, df.values, fmt=fmt)

    # Fill the input namelist
    dest = os.path.join(
        run_info['SIM_INPUT'],
        'input.namelist'
    )
    with open(run_info['INPUT_NAMELIST']) as f:
        newText = f.read().replace('NDAYS', str(int(nr_days)))

    with open(dest, "w") as f:
        f.write(newText)

    # Manipulate some files depending on INPUT_INTERVAL
    stat2bin_src = f"{run_info['TERRA_CFG']}/stat2bin_{run_info['INPUT_INTERVAL']}.f90"
    stat2bin_dst = f"{run_info['TERRA_CFG']}/stat2bin.f90"
    shutil.copy(stat2bin_src, stat2bin_dst)

    # Compile the .dat files into .bin
    for filename in [f"{run_info['TERRA_CFG']}/Makefile",
                     f"{run_info['TERRA_CFG']}/stat2bin.f90"]:
        shutil.copy(filename,run_info['SIM_INPUT'])

    # Compile
    os.chdir(run_info['SIM_INPUT'])
    os.system('make clean &&'
              'make &&'
              './stat2bin')

    return


def create_params(run_info):

    '''
    adds paramater file

    Parameters
    ----------
    fname (string)       : name of forcing file
    info (dictionary)    : information from forcing file
    sitedata (dataframe) : site data table
    '''

    # Read the data
    ds = xr.open_dataset(run_info['MET_FORCING'])

    # Get number of days according to TERRA (longer than forcing).
    # Add one minute, to end date, to make sure that end dates ending on
    # 00:00:00 are considered as a new day.
    nr_days = abs(
        (pd.to_datetime(run_info['START_DATE']) -
         (pd.to_datetime(ds.attrs['time_coverage_end']) + timedelta(minutes=1))).days
    )
    print(nr_days)

    # Set number of required timesteps by TERRA: days * hours * model timestep
    timesteps = (nr_days * 24 * 60)
    time_coverage_end = ds.attrs['time_coverage_end']
    time_start = run_info['START_DATE']

    year       = str(int(time_start.strftime('%Y')))
    month      = str(int(time_start.strftime('%m'))).zfill(2)
    day        = str(int(time_start.strftime('%d'))).zfill(2)
    hour       = str(int(time_start.strftime('%H'))).zfill(2)
    RTydate_ini= f'{year}{month}{day}{hour}'

    # Append info to sitedatafile, store as new params file
    sim_info = pd.DataFrame(columns=['value'])

    # Internal model timestep fixed at 60 seconds.
    sim_info.loc['RTdt', 'value'] = float(60)
    sim_info.loc['RTydate_ini', 'value'] = RTydate_ini
    sim_info.loc['RTnstep_max', 'value'] = int(timesteps)
    sim_info.loc['NDAYS', 'value'] = nr_days # Store, these are days simulated by terra.
    sim_info.loc['START_RUN', 'value'] = run_info['START_RUN']
    sim_info.loc['START_DATE', 'value'] = time_start  # Start run, should be same as RTydate_ini
    sim_info.loc['timestep_interval_seconds', 'value'] = ds.attrs['timestep_interval_seconds']
    sim_info.loc['time_coverage_start', 'value'] = ds.attrs['time_coverage_start']
    sim_info.loc['time_coverage_end', 'value'] = time_coverage_end
    sim_info.loc['time_analysis_start', 'value'] = ds.attrs['time_analysis_start']
    sim_info.loc['local_utc_offset_hours', 'value'] = ds.attrs['local_utc_offset_hours']

    # Fractions
    imp   = float(ds['impervious_area_fraction'])
    tree  = float(ds['tree_area_fraction'])
    grass = float(ds['grass_area_fraction'])
    bare  = float(ds['bare_soil_area_fraction'])
    water = float(ds['water_area_fraction'])

    # Reset removing water fraction
    imp   = imp * (1/(1-water))
    tree  = tree * (1/(1-water))
    grass = grass * (1/(1-water))
    bare  = bare * (1/(1-water))

    # Rescale natural fractions to unity.
    veg   = tree + grass + bare
    tree  = tree * (1/veg)
    grass = grass * (1 / veg)
    bare  = bare * (1 / veg)

    sim_info.loc['EPCOV', 'value'] = np.round(tree+grass,3)
    sim_info.loc['urb_weight', 'value'] = imp
    sim_info.loc['veg_weight', 'value'] = veg
    sim_info.loc['tree', 'value'] = tree
    sim_info.loc['grass', 'value'] = grass
    sim_info.loc['bare', 'value'] = bare

    # Add other required parameters embedded within netcdf
    add_from_nc = [
        "topsoil_clay_fraction",
        "topsoil_sand_fraction",
        "building_mean_height",
        "canyon_height_width_ratio",
        "roof_area_fraction",
        "average_albedo_at_midday",
        f"displacement_height{run_info['DISP_ROUGH']}",
        f"roughness_length_momentum{run_info['DISP_ROUGH']}",
        "anthropogenic_heat_flux_mean",
        "latitude",
        "longitude",
        "measurement_height_above_ground",
    ]
    for item in add_from_nc:
        sim_info.loc[item, 'value'] = float(ds[item])

    #if not os.path.exists(run_info['SIM_INFO']):
    sim_info.to_csv(run_info['SIM_INFO'])

    return


def _replace_mutiple(string, substitutions):

    substrings = sorted(substitutions, key=len, reverse=True)
    regex = re.compile('|'.join(map(re.escape, substrings)))
    return regex.sub(lambda match: substitutions[match.group(0)], string)


# Source SURY: https://github.com/hendrikwout/sury/blob/master/sury.py
def _sury(alb=0.101, \
         emi=0.86, \
         lam_s=0.767, \
         Cv_s=1.25e6, \
         h=15., \
         htw=1.5, \
         roof_f=0.667, \
         z=0., \
         lam_soil=0.28, \
         Cv_soil=1.35e6, \
         alb_snow=0.70, \
         emi_snow=0.997, \
         snow_f=0.0, \
         ustar=0.25, \
         **kwargs):
    """
    Title:   SURY: the Semi-empirical URban canopY parametrization
    Purpose: SURY converts urban canopy parameters - containing the three-dimensional
	     information such as from WUDAPT - into bulk parameters. It can be used
	     with existing bulk urban land-surface schemes for providing canopy-
	     dependent urban physics in atmospheric modelling with a low computation cost.

    Reference:
             Wouters, H., Demuzere, M., Blahak, U., Fortuniak, K., Maiheu., B.,
             Camps, J., Tielemans, and N. P. M. van Lipzig, 2016.  The efficient
             urban-canopy dependency parametrization SURY (v1.0) for atmospheric modelling:
             description and application with the COSMO-CLM model (v5.0_clm6) for a
             Belgian Summer, Geosci. Model Dev., 2016.
    Version: 1.0
    Author: Hendrik Wouters <hendrik.wouters@kuleuven.be>

    License: GNU GENERAL PUBLIC LICENSE Version 3 (license is provided along)
    Input (default values are listed above):
        alb:      albedo [ 1 ]
        emi:      emissivity [ 1 ]
        lam_s:    surface heat conductivity [W m-1 K-1]
        Cv_s:     surface heat capacity [J m-3 K-1]
        H:        building height [m]
        htw:      canyon height-to-width ratio [1],  validity range: 0..2
        roof_f:   roof fraction [ 1 ]
        z:        depths of the ground layers in the vertical column of
                  the bulk land-surface module. Please note that this can be specified
                  as a list or array, eg., [0.01,0.035,0.08,0.17,0.35,0.71,1.43,2.87,5.75,11.51]
        lam_soil: heat conductivity of the soil below [W m-1 K-1]
        Cv_soil:  heat capacity of the soil below [J m-3 K-1]
        alb_snow: albedo of snow [ 1 ]
        emi_snow: emissivity of snow [ 1 ]
        snow_f:   snow fraction [ 1 ]
        ustar:    friction velocity [ m s-1 ]
        --- When provided, Cv_s is obtained from the weighted average of the
            values $C_v,i} and $\lambda_i$) for wall, roof and road according to their
            surface fractions in the urban canopy (see Eqns. 10 of Wouters et al., 2016)
        Cv_roof: heat capacity of the roofs [J m-3 K-1]
        Cv_wall: heat capacity of the walls [J m-3 K-1]
        Cv_road: heat capacity of the roads [J m-3 K-1]
        --- When provided, lam_s is obtained from the weighted average of the
            values $C_v,i} and $\lambda_i$) for wall, roof and road according to their
            surface fractions in the urban canopy, see Eq.  11 of Wouters et al., 2016)
        lam_roof: heat conductivity of the roofs [W m-1 K-1]
        lam_wall: heat conductivity of the walls [W m-1 K-1]
        lam_road: heat conductivity of the roads [W m-1 K-1]
        --- When provided, alb is obtained from the
            values  for wall, roof and road according to  Eq. 16 of Wouters et al. (2016)
        alb_roof: heat conductivity of the roofs [W m-1 K-1]
        alb_wall: heat conductivity of the walls [W m-1 K-1]
        alb_road: heat conductivity of the roads [W m-1 K-1]
        --- When provided, emi is obtained from the
            values  for wall, roof and road analogous to  Eq. 16 of Wouters et al. (2016)
        emi_roof: heat conductivity of the roofs [W m-1 K-1]
        emi_wall: heat conductivity of the walls [W m-1 K-1]
        emi_road: heat conductivity of the roads [W m-1 K-1]
    Output:
        alb_bulk: bulk albedo [ 1 ]
        emi_bulk: bulk emissivity [ 1 ]
        lam_bulk: vertical profile bulk heat conductivity [W m-1 K-1]
        Cv_bulk: vertical profile substrate heat capacity [J m-3 K-1]
        z0 : aerodynamic roughness length [ m ]
        kBm1: kB^{-1} = log(z0/z0H), where z0H is the thermal roughness length [ 1 ]
    Usage:
        import sury

        # Use the default urban-canopy parameters:
        sury.sury()
        # Some user-specified urban-canopy parameters (As always in Python, please explicitly
        # specify the parameter names when deviating from the default argument order):
        sury.sury(h = 20., htw = 10., roof_f = 0.5)
    """
    from numpy import exp, amin, zeros_like

    if htw > 2.0:
        print(
            'Warning, specified canyon height-to-width-ratio falls outside of validity range. In order to extend the formulation, please contact the author.')

    # surface area index of a perfect canyon: parallel canyons and flat roofs
    SAI = (1. + 2. * htw) * (1. - roof_f) + roof_f

    if ('Cv_roof' in kwargs.keys()) or ('Cv_wall' in kwargs.keys()) or ('Cv_road' in kwargs.keys()):
        Cv_s = (kwargs['Cv_roof'] * roof_f + (kwargs['Cv_wall'] * 2. * htw + kwargs['Cv_road']) * (1. - roof_f)) / (
                    roof_f + (2. * htw + 1.) * (1. - roof_f))
        print('Warning, Cv_s is calculated from the values for heat capacity from roads, walls and roofs.' + str(Cv_s))
    if ('lam_roof' in kwargs.keys()) or ('lam_wall' in kwargs.keys()) or ('lam_road' in kwargs.keys()):
        lam_s = (kwargs['lam_roof'] * roof_f + (kwargs['lam_wall'] * 2. * htw + kwargs['lam_road']) * (1. - roof_f)) / (
                    roof_f + (2. * htw + 1.) * (1. - roof_f))
        print(
            'Warning, lam_s is calculated from the values for heat capacity from roads, walls and roofs: ' + str(lam_s))

    # surface-level bulk thermal properties
    Cv_bulk_s = SAI * Cv_s
    lam_bulk_s = SAI * lam_s

    # vertical profiles for bulk heat conductivity and heat capacity

    weight = amin([z, zeros_like(z) + h], axis=0)
    Cv_bulk = (1. - weight / h) * Cv_bulk_s + weight / h * Cv_soil
    lam_bulk = (1. - weight / h) * lam_bulk_s + weight / h * lam_soil

    # bulk radiative properties
    psi_canyon = exp(-0.6 * htw)
    psi_bulk = roof_f + (1. - roof_f) * psi_canyon

    if ('alb_roof' in kwargs.keys()) or ('alb_wall' in kwargs.keys()) or ('alb_road' in kwargs.keys()):
        alb_roof_snow = kwargs['alb_roof'] * (1. - snow_f) + alb_snow * snow_f
        alb_road_snow = kwargs['alb_road'] * (1. - snow_f) + alb_snow * snow_f
        alb_wall_snow = kwargs['alb_wall'] * (1. - snow_f) + alb_snow * snow_f
        alb_bulk = (alb_road_snow + 2. * htw * alb_wall_snow) / (1. + 2. * htw) * psi_canyon * (
                    1. - roof_f) + alb_roof_snow * roof_f

    else:
        alb_bulk = ((1. - snow_f) * alb + snow_f * alb_snow) * psi_bulk

    if ('emi_roof' in kwargs.keys()) or ('emi_wall' in kwargs.keys()) or ('emi_road' in kwargs.keys()):
        albth_roof_snow = (1. - kwargs['emi_roof']) * (1. - snow_f) + (1. - alb_snow) * snow_f
        albth_road_snow = (1. - kwargs['emi_road']) * (1. - snow_f) + (1. - emi_snow) * snow_f
        albth_wall_snow = (1. - kwargs['emi_wall']) * (1. - snow_f) + (1. - emi_snow) * snow_f
        emi_bulk = 1. - ((albth_road_snow + 2. * htw * albth_wall_snow) / (1. + 2. * htw) * psi_canyon * (
                    1. - roof_f) + albth_roof_snow * roof_f)

    else:
        emi_bulk = 1. - psi_bulk * (1. - ((1. - snow_f) * emi + snow_f * emi_snow))

    # surface-layer turbulence properties
    z_0 = 0.075 * h
    nu = 1.461e-5
    Re = ustar * z_0 / nu
    kBm1 = 1.29 * Re ** 0.25 - 2.

    return {'alb_bulk': alb_bulk, \
            'emi_bulk': emi_bulk, \
            'lam_bulk': lam_bulk, \
            'Cv_bulk': Cv_bulk, \
            'z_0': z_0, \
            'kBm1': kBm1 \
            }


def _get_soil_type_class(run_info) -> int:

    sim_info = pd.read_csv(run_info['SIM_INFO'], index_col=0, sep=',')
    topsoil_clay_fraction = float(sim_info.loc['topsoil_clay_fraction','value'])
    topsoil_sand_fraction = float(sim_info.loc['topsoil_sand_fraction','value'])

    # Read TERRA look-up table for soil types.
    SOIL_DATA = pd.read_csv(f"{run_info['LOOKUPDIR']}/data_soil.csv", index_col=0)

    # Do not consider the following types: ICE, ROCK, sea water, sea ice, urban
    SOIL_DATA = SOIL_DATA.loc[
                :,
                ['sand', 'sandy loam', 'loam', 'clay loam', 'clay', 'peat']]

    # Get sand and clay fractions
    csandf = SOIL_DATA.loc['csandf',:]/100 # Re-scale to match site-data
    cclayf = SOIL_DATA.loc['cclayf',:]/100 # Re-scale to match site-data

    # Get differences between site data and look-up table
    diff_csandf = abs(topsoil_sand_fraction-csandf)
    diff_cclayf = abs(topsoil_clay_fraction-cclayf)

    # Find the soiltype index with lowers differences for sand and clay combined
    diff_sum = diff_csandf+diff_cclayf

    # Get index/indices of minimum values
    # If more than 1 index, just pick the first one.
    index_min_diff_sum = np.where(diff_sum == diff_sum.min())[0].tolist()[0]

    # Return index, yet add 1 (.f90 1-based) and 2 as first two columns were skipped
    return int(index_min_diff_sum) + 1 + 2


def put_input_params(run_info):

    sim_info = pd.read_csv(run_info['SIM_INFO'],index_col=0,sep=',')

    # Water puddle storage values from Wouters et al., 2015
    DELTA_M  = 0.12 # delta_m = 12% = 0.12
    WIMPC_C  = np.round(DELTA_M**(-3./2.),2)
    WIMPM_C  = 1.31/1e3   # the maximum water storage w_m = 1.31 [kg/m2]

    # Emis and Heat Cond/Cap from Table 1 (bulk) Wouters et al., 2016
    CSALBVAL = 0.081
    CTALBVAL = 0.89
    CALAVAL  = 1.55
    CRHOCVAL = 2.5E6

    # Approximate LAI min/max by averaging EXTPAR LAI values for globcover's
    # tree (classes 5-10), grass (classes 11-12) and bare soil (class 20), and
    # multiply those with the unity fractions for tree, grass and bare provided
    # in the metadata.
    GC_DATA = pd.read_csv(f"{run_info['LOOKUPDIR']}/globcover_lookup.csv", index_col=1)
    LAI_MX_TREE = GC_DATA.loc[5:10,'LAI_MX'].mean()
    LAI_MX_GRASS = GC_DATA.loc[11:12, 'LAI_MX'].mean()
    LAI_MX_BARE = GC_DATA.loc[20, 'LAI_MX'].mean()
    LAI_MN_TREE = GC_DATA.loc[5:10,'LAI_MN'].mean()
    LAI_MN_GRASS = GC_DATA.loc[11:12, 'LAI_MN'].mean()
    LAI_MN_BARE = GC_DATA.loc[20, 'LAI_MN'].mean()

    LAI_MX = float(sim_info.loc['tree','value']) * LAI_MX_TREE + \
             float(sim_info.loc['grass','value']) * LAI_MX_GRASS + \
             float(sim_info.loc['bare','value']) * LAI_MX_BARE

    LAI_MN = float(sim_info.loc['tree','value']) * LAI_MN_TREE + \
             float(sim_info.loc['grass','value']) * LAI_MN_GRASS + \
             float(sim_info.loc['bare','value']) * LAI_MN_BARE
    ROOTDP = 2.5

    # Set other parameters different for urb and veg ...
    if run_info['CASE'] == "urb":

        EPCOVMA = EPCOVMI = 0
        EZOC = 0.8
        EST = 11
        ELAIMA = ELAIMI = 0
        EKBMTC = 1
        # As used in Demuzere et al. (2017)
        SITSOIL0 = SITSOIL1 = SITSOIL2 = SITSOIL3 = SITSOIL4 = 288.8
        SIWG0 = SIWG1 = SIWG2 = 0.002

    else:
        EPCOVMA = EPCOVMI = pd.to_numeric(sim_info.loc['EPCOV','value'])
        EZOC = 0.08
        EST = _get_soil_type_class(run_info)
        ELAIMA = np.round(LAI_MX,2)
        ELAIMI = np.round(LAI_MN,2)
        EKBMTC = 0
        # Values from Lindenberg setup
        SITSOIL0 = 274.0
        SITSOIL1 = 275.0
        SITSOIL2 = 277.0
        SITSOIL3 = 279.0
        SITSOIL4 = 280.0

        # For initial soil moisture, lower values for dry sites using SIGSCALER
        site_lookup = pd.read_csv(os.path.join(
            f"{run_info['LOOKUPDIR']}",
            "site_lookup.csv"),
            index_col=0)
        SIGW_SCALER = int(site_lookup.loc[run_info['SITE'],'SIGW_SCALER'])
        print(f"> Initial soil moisture scaled by {SIGW_SCALER}")
        SIWG0 = 0.225 / SIGW_SCALER
        SIWG1 = 0.243 / SIGW_SCALER
        SIWG2 = 0.226 / SIGW_SCALER


    # Set type of output here.
    if run_info['OUTPUT_TYPE'] == 'GRB': #OUTPUT IN GRIB
        NTYPE_OUTPUT = 3
        OUTPREFIX = ''
    elif run_info['OUTPUT_TYPE'] == 'ASC': #OUTPUT IN ASCII
        NTYPE_OUTPUT = 5
        OUTPREFIX = 'output.csv'

    # For detailed, take albedo and roughness from sim_info
    # Use H, roof_f and htw as input to SURY.
    if run_info['EXP'] == "d":

        print("=========================================")
        print("EXP = d, using DETAILED information")
        print("=========================================")

        # Results from SURY and directly from sim_info
        sury_out = _sury(
            H=pd.to_numeric(sim_info.loc['building_mean_height','value']),
            htw=pd.to_numeric(sim_info.loc['canyon_height_width_ratio','value']),
            roof_f=pd.to_numeric(sim_info.loc['roof_area_fraction','value'])
        )

        CTALBVAL = 1-sury_out['emi_bulk']
        CRHOCVAL = sury_out['Cv_bulk']
        CALAVAL = sury_out['lam_bulk']
        # CSALBVAL = sury_out['alb_bulk']
        # EZOC = sury_out['z_0']
        # print(f"---->>> EZOC from SURY: {EZOC}")
        CSALBVAL = pd.to_numeric(sim_info.loc['average_albedo_at_midday','value'])
        EZOC  = pd.to_numeric(
            sim_info.loc[
                f"roughness_length_momentum{run_info['DISP_ROUGH']}",
                'value']
        )

    # When lcz is selected, take all params from Stewart et al look-up.
    # See also: https://github.com/matthiasdemuzere/WUDAPT-to-COSMO
    # Surface fractions from LCZ are not explicitly used.
    elif run_info['EXP'] == "lcz" and run_info['CASE'] == "urb":

        print("=========================================")
        print("EXP = lcz, using LCZ-BASED information")
        print("=========================================")

        ucp = pd.read_csv(run_info['UCP_FILE'],
                          sep=';', index_col=0).iloc[:17, :]

        # Read urban canopy parameters ...
        data = ucp.loc[int(pd.to_numeric(run_info['LCZ']))].to_dict()

        # Get the results ...
        sury_out = _sury(
            H=data['URB_BLDH'],
            htw=data['URB_H2W'],
            roof_f=data['URB_BLDFR'],
            Cv_roof=data['URB_RfHCAP'],
            Cv_wall=data['URB_WaHCAP'],
            Cv_road=data['URB_RdHCAP'],
            lam_roof=data['URB_RfHCON'],
            lam_wall=data['URB_WaHCON'],
            lam_road=data['URB_RdHCON'],
            alb_roof=data['URB_RfALB'],
            alb_wall=data['URB_WaALB'],
            alb_road=data['URB_RdALB'],
            emi_roof=data['URB_RfEMI'],
            emi_wall=data['URB_WaEMI'],
            emi_road=data['URB_RdEMI']
        )

        CTALBVAL = 1-sury_out['emi_bulk']
        CSALBVAL = sury_out['alb_bulk']
        CRHOCVAL = sury_out['Cv_bulk']
        CALAVAL = sury_out['lam_bulk']
        EZOC = sury_out['z_0']

    input_array_dict = \
        {# RUN_TERRA
         'RTdt': int(pd.to_numeric(sim_info.loc['RTdt','value'])),
         'RTydate_ini': sim_info.loc['RTydate_ini','value'],
         'RTnstep_max': int(pd.to_numeric(sim_info.loc['RTnstep_max','value'])),
         # EXTPARA
         'EST': EST,
         'EPCOVMA': EPCOVMA,
         'EPCOVMI': EPCOVMI,
         'ROOTDP': ROOTDP,
         'ELAIMA': ELAIMA,
         'ELAIMI': ELAIMI,
         'EZ0C': EZOC,
         'EKBMTC': EKBMTC,
         'ELAT': pd.to_numeric(sim_info.loc['latitude','value']),
         'WIMPM_C': WIMPM_C,
         'WIMPC_C': WIMPC_C,
        # SOILINIT
         'SITSOIL0': SITSOIL0,
         'SITSOIL1': SITSOIL1,
         'SITSOIL2': SITSOIL2,
         'SITSOIL3': SITSOIL3,
         'SITSOIL4': SITSOIL4,
         'SIWG0': SIWG0,
         'SIWG1': SIWG1,
         'SIWG2': SIWG2,
        # METFORCING
         'MFDZ': pd.to_numeric(sim_info.loc['measurement_height_above_ground','value']),
         'METDIR': run_info['SIM_INPUT'],
        # OUTPUT
         'NTYPE_OUTPUT': int(NTYPE_OUTPUT),
         'OUTDIR': run_info['OUTPUTDIR'],
         'OUTPREFIX': OUTPREFIX,
         'ONOUTINT': run_info['OUTPUT_INTERVAL'],
        }

    data_soil_dict = \
        {'CRHOCVAL': CRHOCVAL,
         'CALAVAL': CALAVAL,
         'CTALBVAL': CTALBVAL,
         'CSALBVAL': CSALBVAL}

    # Convert values to strings
    input_array_dict = {key: str(value) for key, value in input_array_dict.items()}
    data_soil_dict = {key: str(value) for key, value in data_soil_dict.items()}

    # Replace for INPUT_ARRAY
    f1 = str(open(run_info['INPUT_ARRAY']).read());
    input_array_dest = os.path.join(run_info['SIM_CFG'], 'INPUT_TERRA')
    fnew = _replace_mutiple(f1,input_array_dict);
    with open(input_array_dest, "w") as f2:
        f2.write(fnew)

    # Replace for data_soil.f90
    f1 = str(open(run_info['INPUT_DATA_SOIL']).read());
    data_soil_dest = os.path.join(run_info['SRCDIR'],'data_soil.f90')
    fnew = _replace_mutiple(f1,data_soil_dict);
    with open(data_soil_dest, "w") as f_f90:
        f_f90.write(fnew)

    # Also keep track of this file in config
    data_soil_arc = os.path.join(run_info['SIM_CFG'],'data_soil.f90')
    with open(data_soil_arc, "w") as f_arc:
        f_arc.write(fnew)

    if not run_info['EXP'] in ['b','d','lcz']:
        shutil.copy(run_info['INPUT_ARRAY'].replace(".temp",f"_{run_info['CASE']}_{run_info['EXP']}"),
                    os.path.join(run_info['SIM_CFG'], 'INPUT_TERRA'))
        shutil.copy(run_info['INPUT_DATA_SOIL'].replace(".temp", f"_{run_info['EXP']}"),
                    os.path.join(run_info['SRCDIR'],'data_soil.f90'))


def compile_model(run_info):

    print("COMPILING TERRA TSA ...")
    # Choose appropriate terra_io depending on 30 or 60 minute forcing
    terra_io_src = f"{run_info['SRCDIR']}/terra_io_{run_info['INPUT_INTERVAL']}.f90"
    terra_io_dst = f"{run_info['SRCDIR']}/terra_io.f90"
    shutil.copy(terra_io_src, terra_io_dst)

    # Re-compile terra
    os.chdir(run_info['SRCDIR'])
    comp_out = os.system('make clean && make')
    if comp_out !=0:
        print(">>> Problem with compilation ... Abort")
        sys.exit()
    else:
        print(">>> Compilation completed without error.")

    exec_terra_src = os.path.join(run_info['SRCDIR'],"terra")
    dest_terra_src = os.path.join(
        run_info['SIM_CFG'],
        f"terra_{run_info['SITE']}_{run_info['CASE']}_{run_info['EXP']}"
    )
    shutil.copy(exec_terra_src,dest_terra_src)

    print(f"Executable available: {dest_terra_src}")


def run_model(run_info):

    start_time = time.time()

    print("****************************")
    print(f'TERRA simulation started, be patient ...')
    print("****************************")

    os.chdir(run_info['SIM_CFG'])
    os.system(f"./terra_{run_info['SITE']}_{run_info['CASE']}_{run_info['EXP']}")

    end_time = time.time()
    run_time = np.round((end_time - start_time),2) # Time in seconds.

    print("****************************")
    print(f'Executed in {run_time} seconds')
    print("****************************")

    # Write simulation time to sim_info.
    sim_info = pd.read_csv(run_info['SIM_INFO'],index_col=0)
    sim_info.loc[f"RUNTIME_SEC_{run_info['CASE']}", 'value'] = run_time
    sim_info.to_csv(run_info['SIM_INFO'])

    # Remove YUDRDAT file
    os.remove(run_info['YUDRDAT'])


def convert_to_netcdf(run_info):

    sim_info = pd.read_csv(run_info['SIM_INFO'], index_col=0, sep=',')

    ini_date = pd.to_datetime(sim_info.loc['START_DATE', 'value'])
    end_date = pd.to_datetime(sim_info.loc['time_coverage_end', 'value'])

    # Get the months to convert
    months = pd.period_range(ini_date, end_date, freq='M').strftime('%Y%m')

    os.chdir(run_info['OUTPUTDIR'])

    for yearmm in months:

        os.system(f"cat {yearmm}?????? > {yearmm}.grb")
        os.system(f"ncl_convert2nc {yearmm}.grb")

    final_file = f"{run_info['OUTPUTDIR']}/" \
                 f"{run_info['SITE']}_{run_info['CASE']}_{run_info['EXP']}.nc"

    for yearmm in months:

        os.system(f"cdo -f nc cat {yearmm}.nc {final_file}")

    # Remove all but final and info file
    for file in glob.glob(f"{run_info['OUTPUTDIR']}/*"):
        if not file in [final_file,run_info['SIM_INFO']]:
            os.remove(file)

    return final_file


def _geo_idx(dd, dd_array):

   geo_idx = (np.abs(dd_array - dd)).argmin()
   return geo_idx

def _get_ahf_point(ahf_file, lat, lon):

    nc_ahf = nc.Dataset(ahf_file)
    lats = nc_ahf.variables['lat'][:]
    lons = nc_ahf.variables['lon'][:]

    lat_idx = _geo_idx(float(lat), lats)
    lon_idx = _geo_idx(float(lon), lons)

    ahf = nc_ahf['AHF'][lat_idx, lon_idx].data.max()

    return ahf

# Eq. A7 Flanner et al. (2009)
def _get_ahf_amplitude(lat):

    if lat > 33:
        A2 = 1 - np.e ** -((lat-33)/25)
    elif lat < -33:
        A2 = -(1 - np.e ** -((lat + 33) / 25))
    else:
        A2 = 0

    return A2

def _get_ahf_timeseries(ahf, lat, td, ty):

    """
    mean annual AHF values are provided, and scaled with weighting functions
    dependent on local time of day and time of year.

    Original paper (Flanner, 2009) contains some errata. These are fixed in this version:
    http://clasp-research.engin.umich.edu/faculty/flanner/content/ppr/ppr_ener.pdf
    """

    # Parameters from Table 3 in Flanner (2009)
    ahf_params = {
        'b1': 0.451,     # Daytime amplification
        'b2': 0.8,       # Daytime amplification
        'sigma': 0.18,   # Slope of evening tapering
        'mu': 0.5,       # Timing of morning and evening ramps
        'A1': -0.3,      # Amplitude of morning and evening traffic maxima
        'f': 2.0,        # Timing of morning and evening traffic maxima
        'alpha': 10.0,   # Slope of morning increase
        'epsilon': 0.25, # Timing of morning and evening ramps
    }

    # Eqs. A4/5
    erf_e1_arg = ahf_params['alpha'] * \
                 (td - ahf_params['mu'] + ahf_params['epsilon']) / \
                 ahf_params['sigma']
    E1 = 0.5 * (erf(erf_e1_arg) + 1.0)

    erf_e2_arg = -ahf_params['alpha'] * \
                 (td - ahf_params['mu'] - ahf_params['epsilon']) / \
                 ahf_params['sigma']
    E2 = 0.5 * (erf(erf_e2_arg) + 1.0)

    # Eq. A3
    H = ahf_params['A1'] * np.cos(2 * np.pi * ahf_params['f'] * td)

    # Eq. A2
    N = 1 / (ahf_params['sigma'] * np.sqrt(2*np.pi)) * \
        np.e**( -(td - ahf_params['mu'])**2 / (2*ahf_params['sigma']**2) )

    # Eq. A1
    wd = (N * E1 * ahf_params['b1']) + (H * E1 * E2) + ahf_params['b2']
    #plt.plot(td,wd)

    # Eq. A7
    A2 = _get_ahf_amplitude(float(lat))

    # Eq. A6
    wy = 1 + A2 * np.sin( 2 * np.pi * (ty + 0.25))
    #plt.plot(wy)

    # Get timeseries
    ahf_series = ahf * wy * wd

    return ahf_series


def create_final_df_from_grb(run_info, sim_info):

    '''
    TERRA input is always extended to last (half) hour of the day.
    First make dataframe to that time step, then crop using 'time_coverage_end'
    '''

    #ini_date = pd.to_datetime(sim_info['time_coverage_start'])
    ini_date = pd.to_datetime(sim_info['START_DATE'])
    end_date = pd.to_datetime(sim_info['time_coverage_end'])

    timestep_interval = sim_info['timestep_interval_seconds']
    #print(ini_date, end_date, timestep_interval)

    # Get the dates as used in the simulation
    dates = pd.date_range(
        ini_date,
        end_date,
        freq=f"{pd.to_numeric(timestep_interval)/3600}H"
    )

    # Get the dates used in the simulation (
    # Get the forcing variables
    # load forcing information from netcdf
    ds = xr.open_dataset(run_info['MET_FORCING'])

    # Get local time offset, needed for the AHF cycle dates
    local_utc_offset_hours = ds.attrs['local_utc_offset_hours']

    # create empty dataframe, to into account debugging.
    forcing = pd.DataFrame(
        index=dates,
        columns=['SWdown', 'LWdown', 'Tair', 'Qair', 'PSurf', 'Wind'])

    for var_i, var in enumerate(forcing.columns.to_list()):

        if var != 'Wind':
            forcing.loc[:,forcing.columns[var_i]] = \
                ds[var].sel(time=slice(ini_date,end_date)).values.flatten()
        else:
            Wind_E = ds['Wind_E'].sel(time=slice(ini_date,end_date)).values.flatten()
            Wind_N = ds['Wind_N'].sel(time=slice(ini_date, end_date)).values.flatten()
            Wind = np.sqrt(Wind_E ** 2 + Wind_N ** 2)
            forcing.loc[:, 'Wind'] = Wind

    # Read the raw results
    OUTPUT_BASE = run_info['OUTPUTDIR'][:-3] # drop urb or veg
    urb_nc = xr.open_dataset(f"{OUTPUT_BASE}/urb/"
                             f"{run_info['SITE']}_urb_{run_info['EXP']}.nc")
    veg_nc = xr.open_dataset(f"{OUTPUT_BASE}/veg/"
                             f"{run_info['SITE']}_veg_{run_info['EXP']}.nc")
    # Get the urb and veg fractions
    urb_w = float(sim_info['urb_weight'])
    veg_w = float(sim_info['veg_weight'])

    # Get the mean annual anthropogenic heat flux
    if run_info['EXP'] == 'b':
        ahf = _get_ahf_point(
            run_info['AHF_FILE'],
            sim_info['latitude'],
            sim_info['longitude']
        )
    elif run_info['EXP'] == 'd':
        ahf = float(sim_info['anthropogenic_heat_flux_mean'])
    elif run_info['EXP'] == 'lcz':
        ucp = pd.read_csv(
            run_info['UCP_FILE'],
            sep=';',
            index_col=0).iloc[:17, :]
        ahf = ucp.loc[int(pd.to_numeric(run_info['LCZ']))].to_dict()['AHF']

    # Get fractional hours and days of year
    offset    = pd.Timedelta('%s hours' % local_utc_offset_hours)
    dates_ahf = dates + offset

    if int(run_info['INPUT_INTERVAL']) == 30:
        td = (dates_ahf.hour + dates_ahf.minute/60)/24
    elif int(run_info['INPUT_INTERVAL']) == 60:
        td = dates_ahf.hour/24
    ty = dates_ahf.dayofyear/366

    # Get the AHF time series corresponding to model output.
    ahf_series = _get_ahf_timeseries(ahf, sim_info['latitude'], td, ty)

    # variables of interest, and available from TERRA_URB
    # Units as comments: Urban Plumber | TERRA out
    vars_dict_to_weight = \
        {'SWnet'   : 'ASOB_S_GDS10_SFC',      # W/m2 | W/m2
         'LWnet'   : 'ATHB_S_GDS10_SFC',      # W/m2 | W/m2
         'Qle'     : 'ALHFL_S_GDS10_SFC',     # W/m2 | W/m2
         'Qh'      : 'ASHFL_S_GDS10_SFC',     # W/m2 | W/m2
         'Rainf'   : 'RAIN_GSP_GDS10_SFC',    # kg/m2/s | kg/m2
         'Snowf'   : 'SNOW_GSP_GDS10_SFC',    # kg/m2/s | kg/m2
         }

    # Create dataframe to store all results
    df_column_names = ['SWnet', 'LWnet', 'SWup', 'LWup',
                       'Qle', 'Qh', 'Qanth','Qstor',
                       'Albedo', 'AvgSurfT',
                       'Rainf', 'Snowf', 'Qs', 'Qsb']

    # URB simulation
    df_urb = pd.DataFrame(index=dates,
                          columns=df_column_names)

    for up_key, terra_key in vars_dict_to_weight.items():
        df_urb[up_key] = urb_nc[terra_key].values.flatten()[:len(dates)]

    df_urb['Qanth']    = ahf_series
    df_urb['Qh']       = df_urb['Qh'] + df_urb['Qanth']
    df_urb['Qstor']    = df_urb['SWnet'] + df_urb['LWnet'] + df_urb['Qanth'] \
                                  - df_urb['Qle'] - df_urb['Qh']
    df_urb['SWup']     = forcing['SWdown'] - df_urb['SWnet']
    df_urb['LWup']     = forcing['LWdown'] - df_urb['LWnet']
    df_urb['Albedo']   = df_urb['SWup'] / forcing['SWdown']

    urb_tmp_qs         = urb_nc['RUNOFF_GDS10_DBLY'].values
    df_urb['Qs']       = urb_tmp_qs[:len(dates), 0, 0, 0] #TODO: TERRA in kg/m2, UP expects kg/m2/s?
    df_urb['Qsb']      = urb_tmp_qs[:len(dates), 1, 0, 0] #TODO: TERRA in kg/m2, UP expects kg/m2/s?

    # First level of soil temparture is surface
    urb_t_so           = urb_nc['SO_TEMP_GDS10_DBLL'].values
    df_urb['AvgSurfT'] = urb_t_so[:len(dates), 0, 0, 0]

    # VEG simulation
    df_veg = pd.DataFrame(index=dates,
                          columns=df_column_names)
    for up_key, terra_key in vars_dict_to_weight.items():
        df_veg[up_key] = veg_nc[terra_key].values.flatten()[:len(dates)]

    df_veg['Qanth']    = ahf_series
    df_veg['Qh']       = df_veg['Qh'] + df_veg['Qanth']
    df_veg['Qstor']    = df_veg['SWnet'] + df_veg['LWnet'] + df_veg['Qanth'] \
                                    - df_veg['Qle'] - df_veg['Qh']
    df_veg['SWup']     = forcing['SWdown'] - df_veg['SWnet']
    df_veg['LWup']     = forcing['LWdown'] - df_veg['LWnet']
    df_veg['Albedo']   = df_veg['SWup'] / forcing['SWdown']

    veg_tmp_qs         = veg_nc['RUNOFF_GDS10_DBLY'].values
    df_veg['Qs']       = veg_tmp_qs[:len(dates), 0, 0, 0]
    df_veg['Qsb']      = veg_tmp_qs[:len(dates), 1, 0, 0]

    # First level of soil temparture is surface
    veg_t_so           = veg_nc['SO_TEMP_GDS10_DBLL'].values
    df_veg['AvgSurfT'] = veg_t_so[:len(dates), 0, 0, 0]

    # Put results together with weights
    df = (df_urb * urb_w) + (df_veg * veg_w)

    # Fix some units, in line with UP expectations
    if int(run_info['INPUT_INTERVAL']) == 30:
        df['Rainf'] = df['Rainf'] / 1800
        df['Snowf'] = df['Snowf'] / 1800
        df['Qs'] = df['Qs'] / 1800
        df['Qsb'] = df['Qsb'] / 1800
    elif int(run_info['INPUT_INTERVAL']) == 60:
        df['Rainf'] = df['Rainf'] / 3600
        df['Snowf'] = df['Snowf'] / 3600
        df['Qs'] = df['Qs'] / 3600
        df['Qsb'] = df['Qsb'] / 3600

    # Additional data, that does not need weighing
    #df['Qanth_Qh'] = df['Qanth']

    # Strange values for SnowT??
    #df['SnowT'] = veg_nc.variables['T_SNOW_GDS10_SFC'][:].flatten().data[:len(dates)] # From VEG only.

    # Add to results df
    data = pd.concat([df, forcing], axis=1)

    # Set infinity to NAN.
    data.replace([np.inf, -np.inf], np.nan, inplace=True)

    # Just extract SoilMoist and SoilTemp from netcdf.
    data_SoilTemp  = veg_t_so[:len(dates), 1:, 0, 0]       # K

    veg_w_so = veg_nc['SO_WA_ICE_GDS10_DBLL'].values      # kg/m2
    data_SoilMoist = veg_w_so[:len(dates), :, 0, 0]

    # Save as file - Taking a lot of storage, ignore by default
    #data.to_csv(run_info['FINAL_FNAME_OUT'].replace('.nc','.csv'))

    # In case of a production run, expand time to include all timesteps,
    # starting at 'time_coverage_start'
    if 'production' in sim_info['START_RUN']:

        print(' >> This is a PRODUCTION run, pre-pending NaNs to cover all time steps')

        start_coverage_date = pd.to_datetime(sim_info['time_coverage_start'])
        all_dates = pd.date_range(
            start_coverage_date,
            end_date,
            freq=f"{pd.to_numeric(timestep_interval) / 3600}H"
        )

        # Pre-pend NaNs to beginning of data df
        data = data.reindex(all_dates, fill_value=np.nan)

        # Get the number of NaNs introduced
        nr_nans = np.sum(np.isnan(data.SWnet))

        # Append these NaNs also to _SoilMoist and _SoilTemp
        nan_arr = np.full([nr_nans, 8], np.nan)
        data_SoilMoist = np.append(nan_arr, data_SoilMoist, axis=0)
        data_SoilTemp = np.append(nan_arr, data_SoilTemp, axis=0)

    return data, data_SoilMoist, data_SoilTemp


def create_empty_netcdf(run_info: dict,
                        sim_info: dict,
                        data,
                        num_soil_layers: int,
                        missing_float: float) -> None:

    '''creates empty netcdf dataset complying with Urban-PLUMBER protocol v1.0

    Inputs
    ------
    output filename string.

    '''

    # setting coordinate values
    if 'production' in sim_info['START_RUN']:
        start_date = sim_info['time_coverage_start']
    else:
        start_date = sim_info['START_DATE']

    timesteps = data.shape[0]
    times = [t * int(pd.to_numeric(run_info['OUTPUT_INTERVAL']))*60 for t in range(0, int(timesteps))]
    soil_layers = [i for i in range(1, num_soil_layers + 1)]

    # open netcdf files (r = read only, w = write new)
    with nc.Dataset(filename=run_info['FINAL_FNAME_OUT'], mode='w', format='NETCDF4') as o:

        ############ create dimensions ############
        o.createDimension(dimname='time', size=timesteps)
        o.createDimension(dimname='soil_layer', size=num_soil_layers)
        o.createDimension(dimname='x', size=1)
        o.createDimension(dimname='y', size=1)

        ############ create coordinates ############
        var = 'time'
        o.createVariable(var, datatype='i4', dimensions=('time'), fill_value = missing_float)
        o.variables[var].long_name     = 'Time'
        o.variables[var].standard_name = 'time'
        o.variables[var].units         = 'seconds since %s' %start_date
        o.variables[var].calendar      = 'standard'
        o.variables[var][:]            = times

        var = 'soil_layer'
        o.createVariable(var, datatype='i4', dimensions=('soil_layer'), fill_value = missing_float)
        o.variables[var].long_name     = 'Soil layer number'
        o.variables[var][:]            = soil_layers

        var = 'x'
        o.createVariable(var, datatype='i4', dimensions=('x'), fill_value = missing_float)
        o.variables[var].long_name     = 'x dimension'
        o.variables[var][:]            = 1

        var = 'y'
        o.createVariable(var, datatype='i4', dimensions=('y'), fill_value = missing_float)
        o.variables[var].long_name     = 'y dimension'
        o.variables[var][:]            = 1

        ################### latidude and longitude ###################

        var = 'longitude'
        o.createVariable(var, datatype='f8', dimensions=('y', 'x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Longitude'
        o.variables[var].standard_name = 'longitude'
        o.variables[var].units         = 'degrees_east'
        o.variables[var][:]            = sim_info['longitude']

        var = 'latitude'
        o.createVariable(var, datatype='f8', dimensions=('y', 'x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Latitude'
        o.variables[var].standard_name = 'latitude'
        o.variables[var].units         = 'degrees_north'
        o.variables[var][:]            = sim_info['latitude']

        ##########################################################################
        ################### critical energy balance components ###################

        var = 'SWnet'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Net shortwave radiation (downward)'
        o.variables[var].standard_name = 'surface_net_downward_shortwave_flux'
        o.variables[var].units         = 'W/m2'

        var = 'LWnet'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Net longwave radiation (downward)'
        o.variables[var].standard_name = 'surface_net_downward_longwave_flux'
        o.variables[var].units         = 'W/m2'

        var = 'SWup'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Upwelling shortwave radiation flux (positive upward)'
        o.variables[var].standard_name = 'surface_upwelling_shortwave_flux'
        o.variables[var].units         = 'W/m2'

        var = 'LWup'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Upwelling longwave radiation flux (positive upward)'
        o.variables[var].standard_name = 'surface_upwelling_longwave_flux'
        o.variables[var].units         = 'W/m2'

        var = 'Qle'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Latent heat flux (upward)'
        o.variables[var].standard_name = 'surface_upward_latent_heat_flux'
        o.variables[var].units         = 'W/m2'

        var = 'Qh'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Sensible heat flux (upward)'
        o.variables[var].standard_name = 'surface_upward_sensible_heat_flux'
        o.variables[var].units         = 'W/m2'

        var = 'Qanth'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Anthropogenic heat flux (upward)'
        o.variables[var].standard_name = 'surface_upward_heat_flux_due_to_anthropogenic_energy_consumption'
        o.variables[var].units         = 'W/m2'

        var = 'Qstor'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Net storage heat flux in all materials (increase)'
        o.variables[var].standard_name = 'surface_thermal_storage_heat_flux'
        o.variables[var].units         = 'W/m2'

        ################### additional energy balance compoenents #################

        var = 'Qg'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Ground heat flux (downward)'
        o.variables[var].standard_name = 'downward_heat_flux_at_ground_level_in_soil'
        o.variables[var].units         = 'W/m2'

        var = 'Qanth_Qh'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Anthropogenic sensible heat flux (upward)'
        o.variables[var].standard_name = 'surface_upward_sensible_heat_flux_due_to_anthropogenic_energy_consumption'
        o.variables[var].units         = 'W/m2'

        var = 'Qanth_Qle'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Anthropogenic latent heat flux (upward)'
        o.variables[var].standard_name = 'surface_upward_latent_heat_flux_due_to_anthropogenic_energy_consumption'
        o.variables[var].units         = 'W/m2'

        var = 'Qtau'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Momentum flux (downward)'
        o.variables[var].standard_name = 'magnitude_of_surface_downward_stress'
        o.variables[var].units         = 'N/m2'

        ##################### general water balance components #####################

        var = 'Snowf'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Snowfall rate (downward)'
        o.variables[var].standard_name = 'snowfall_flux'
        o.variables[var].units         = 'kg/m2/s'

        var = 'Rainf'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Rainfall rate (downward)'
        o.variables[var].standard_name = 'rainfall_flux'
        o.variables[var].units         = 'kg/m2/s'

        var = 'Evap'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Total evapotranspiration (upward)'
        o.variables[var].standard_name = 'surface_evapotranspiration'
        o.variables[var].units         = 'kg/m2/s'

        var = 'Qs'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Surface runoff (out of gridcell)'
        o.variables[var].standard_name = 'surface_runoff_flux'
        o.variables[var].units         = 'kg/m2/s'

        var = 'Qsb'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Subsurface runoff (out of gridcell)'
        o.variables[var].standard_name = 'subsurface_runoff_flux'
        o.variables[var].units         = 'kg/m2/s'

        var = 'Qsm'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Snowmelt (solid to liquid)'
        o.variables[var].standard_name = 'surface_snow_and_ice_melt_flux'
        o.variables[var].units         = 'kg/m2/s'

        var = 'Qfz'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Re-freezing of water in the snow (liquid to solid)'
        o.variables[var].standard_name = 'surface_snow_and_ice_refreezing_flux'
        o.variables[var].units         = 'kg/m2/s'

        var = 'DelSoilMoist'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Change in soil moisture (increase)'
        o.variables[var].standard_name = 'change_over_time_in_mass_content_of_water_in_soil'
        o.variables[var].units         = 'kg/m2'

        var = 'DelSWE'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Change in snow water equivalent (increase)'
        o.variables[var].standard_name = 'change_over_time_in_surface_snow_and_ice_amount'
        o.variables[var].units         = 'kg/m2'

        var = 'DelIntercept'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Change in interception storage (increase)'
        o.variables[var].standard_name = 'change_over_time_in_canopy_water_amount'
        o.variables[var].units         = 'kg/m2'

        var = 'Qirrig'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Anthropogenic water flux from irrigation (increase)'
        o.variables[var].standard_name = 'surface_downward_mass_flux_of_water_due_to_irrigation'
        o.variables[var].units         = 'kg/m2/s'

        ########################## surface state variables ########################

        var = 'SnowT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Snow surface temperature'
        o.variables[var].standard_name = 'surface_snow_skin_temperature'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'snow'

        var = 'VegT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Vegetation canopy temperature'
        o.variables[var].standard_name = 'surface_canopy_skin_temperature'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'vegetation'

        var = 'BaresoilT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Temperature of bare soil'
        o.variables[var].standard_name = 'surface_ground_skin_temperature'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'baresoil'

        var = 'AvgSurfT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Average surface temperature (skin)'
        o.variables[var].standard_name = 'surface_temperature'
        o.variables[var].units         = 'K'

        var = 'RadT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Surface radiative temperature'
        o.variables[var].standard_name = 'surface_radiative_temperature'
        o.variables[var].units         = 'K'

        var = 'Albedo'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Surface albedo'
        o.variables[var].standard_name = 'surface_albedo'
        o.variables[var].units         = '1'

        var = 'SWE'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Snow water equivalent'
        o.variables[var].standard_name = 'surface_snow_amount'
        o.variables[var].units         = 'kg/m2'

        var = 'SurfStor'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Surface water storage'
        o.variables[var].standard_name = 'surface_water_amount_assuming_no_snow'
        o.variables[var].units         = 'kg/m2'

        var = 'SnowFrac'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Snow covered fraction'
        o.variables[var].standard_name = 'surface_snow_area_fraction'
        o.variables[var].units         = '1'

        var = 'SAlbedo'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Snow albedo'
        o.variables[var].standard_name = 'snow_and_ice_albedo'
        o.variables[var].units         = '1'
        o.variables[var].subgrid       = 'snow'

        var = 'CAlbedo'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Vegetation canopy albedo'
        o.variables[var].standard_name = 'canopy_albedo'
        o.variables[var].units         = '1'
        o.variables[var].subgrid       = 'vegetation'

        var = 'UAlbedo'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Urban canopy albedo'
        o.variables[var].standard_name = 'urban_albedo'
        o.variables[var].units         = '1'
        o.variables[var].subgrid       = 'urban'

        var = 'LAI'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Leaf area index'
        o.variables[var].standard_name = 'leaf_area_index'
        o.variables[var].units         = '1'
        o.variables[var].subgrid       = 'vegetation'

        var = 'RoofSurfT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Roof surface temperature (skin)'
        o.variables[var].standard_name = 'surface_roof_skin_temperature'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'roof'

        var = 'WallSurfT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Wall surface temperature (skin)'
        o.variables[var].standard_name = 'surface_wall_skin_temperature'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'wall'

        var = 'RoadSurfT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Road surface temperature (skin)'
        o.variables[var].standard_name = 'surface_road_skin_temperature'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'road'

        var = 'TairSurf'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Near surface air temperature (2m)'
        o.variables[var].standard_name = 'air_temperature_near_surface'
        o.variables[var].units         = 'K'

        var = 'TairCanyon'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Air temperature in street canyon (bulk)'
        o.variables[var].standard_name = 'air_temperature_in_street_canyon'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'canyon'

        var = 'TairBuilding'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Air temperature in buildings (bulk)'
        o.variables[var].standard_name = 'air_temperature_in_buildings'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'building'

        ######################## Sub-surface state variables ######################

        var = 'SoilMoist'
        o.createVariable(var, datatype='f8', dimensions=('time','soil_layer','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Average layer soil moisture'
        o.variables[var].standard_name = 'moisture_content_of_soil_layer'
        o.variables[var].units         = 'kg/m2'

        var = 'SoilTemp'
        o.createVariable(var, datatype='f8', dimensions=('time','soil_layer','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Average layer soil temperature'
        o.variables[var].standard_name = 'soil_temperature'
        o.variables[var].units         = 'K'

        ########################## Evaporation components #########################

        var = 'TVeg'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Vegetation transpiration'
        o.variables[var].standard_name = 'transpiration_flux'
        o.variables[var].units         = 'kg/m2/s'
        o.variables[var].subgrid       = 'vegetation'

        var = 'ESoil'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Bare soil evaporation'
        o.variables[var].standard_name = 'liquid_water_evaporation_flux_from_soil'
        o.variables[var].units         = 'kg/m2/s'
        o.variables[var].subgrid       = 'baresoil'

        var = 'RootMoist'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Root zone soil moisture'
        o.variables[var].standard_name = 'mass_content_of_water_in_soil_defined_by_root_depth'
        o.variables[var].units         = 'kg/m2'

        var = 'SoilWet'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Total soil wetness'
        o.variables[var].standard_name = 'relative_soil_moisture_content_above_wilting_point'
        o.variables[var].units         = '1'

        var = 'ACond'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Aerodynamic conductance'
        o.variables[var].standard_name = 'inverse_aerodynamic_resistance'
        o.variables[var].units         = 'm/s'

        ########################## forcing data variables #########################

        var = 'SWdown'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Downward shortwave radiation at measurement height'
        o.variables[var].standard_name = 'surface_downwelling_shortwave_flux_in_air'
        o.variables[var].units         = 'W/m2'

        var = 'LWdown'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Downward longwave radiation at measurement height'
        o.variables[var].standard_name = 'surface_downwelling_longwave_flux_in_air'
        o.variables[var].units         = 'W/m2'

        var = 'Tair'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Air temperature at measurement height'
        o.variables[var].standard_name = 'air_temperature'
        o.variables[var].units         = 'K'

        var = 'Qair'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Specific humidity at measurement height'
        o.variables[var].standard_name = 'surface_specific_humidity'
        o.variables[var].units         = '1'

        var = 'PSurf'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Air pressure at measurement height'
        o.variables[var].standard_name = 'surface_air_pressure'
        o.variables[var].units         = 'Pa'

        var = 'Wind'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Wind speed at measurement height'
        o.variables[var].standard_name = 'wind_speed'
        o.variables[var].units         = 'm/s'

    return


def set_netcdf_data(run_info,
                    data,
                    data_SoilMoist,
                    data_SoilTemp,
                    num_soil_layers,
                    missing_float):

    timesteps = data.shape[0]

    no_data_1D = np.full([timesteps], missing_float)
    no_data_2D = np.full([timesteps, num_soil_layers], missing_float)

    # open netcdf files (r = read only, r+ = append existing)
    with nc.Dataset(filename=run_info['FINAL_FNAME_OUT'], mode='r+', format='NETCDF4') as o:
        # set metadata
        o.title = 'TERRA_URB output for the Urban-PLUMBER project'
        o.site = f"Site name: {run_info['SITE']}"
        o.experiment = f"{run_info['EXP']}"
        o.institution = 'Ruhr University Bochum (Germany) / https://www.climate.ruhr-uni-bochum.de/news/'
        o.primary_contact = 'Matthias Demuzere / Matthias.demuzere@rub.de'
        o.secondary_contact = 'Mikhail Varentsov / mvar91@gmail.com'
        o.model = f"{run_info['TERRA_VERSION']}"
        o.source = 'TERRA_URB, The urban-canopy land-surface scheme of COSMO-CLM'
        o.references = 'http://www.doi.org/10.1016/j.uclim.2014.11.005 ' \
                       'http://www.doi.org/10.5194/gmd-9-3027-2016 ' \
                       'http://www.cosmo-model.org/content/tasks/workGroups/wg3b/docs/terra_urb_user.pdf'
        o.repository = 'https://github.com/matthiasdemuzere/urban-plumber-terra-pub'
        o.site_experience = 'Experience with Toulouse / Melbourne / Singapore / Helsinki: ' \
                            '- Melbourne and Toulouse: http://www.doi.org/10.1002/joc.3656 ' \
                            '- Helsinki: http://www.doi.org/10.1002/qj.2659 ' \
                            '- Singapore: http://www.doi.org/10.1002/qj.3028'
        o.additional_data = 'No additional data has been used'
        o.comment = 'We acknowledge the support of Hendrik Wouters (VITO) and Jan-Peter Schulz (DWD). ' \
                    'Additional information is available in github repo: '\
                    'https://github.com/matthiasdemuzere/urban-plumber-terra-pub.'

        o.history = 'Created with %s at %s' % (__file__, pd.Timestamp.now())

        # Critical energy balance components
        o.variables['SWnet'][:] = data['SWnet'].values  # Net shortwave radiation (downward)
        o.variables['LWnet'][:] = data['LWnet'].values  # Net longwave radiation (downward)
        o.variables['SWup'][:] = data['SWup'].values  # Upwelling shortwave radiation
        o.variables['LWup'][:] = data['LWup'].values  # Upwelling longwave radiation
        o.variables['Qle'][:] = data['Qle'].values  # Latent heat flux (upward)
        o.variables['Qh'][:] = data['Qh'].values  # Sensible heat flux (upward)
        o.variables['Qanth'][:] = data['Qanth'].values  # Anthropogenic heat flux (upward)
        o.variables['Qstor'][:] = data['Qstor'].values  # Net storage heat flux in all materials (increase)
        # Additional energy balance components
        o.variables['Qg'][:] = no_data_1D  # Ground heat flux (downward)
        o.variables['Qanth_Qh'][:] = no_data_1D  # Anthropogenic sensible heat flux (upward)
        o.variables['Qanth_Qle'][:] = no_data_1D  # Anthropogenic latent heat flux (upward)
        o.variables['Qtau'][:] = no_data_1D  # Momentum flux (downward)
        # General water balance components
        o.variables['Snowf'][:] = data['Snowf'].values  # Snowfall rate (downward)
        o.variables['Rainf'][:] = data['Rainf'].values  # Rainfall rate (downward)
        o.variables['Evap'][:] = no_data_1D  # Total evapotranspiration (upward)
        o.variables['Qs'][:] = data['Qs'].values  # Surface runoff (out of gridcell)
        o.variables['Qsb'][:] = data['Qsb'].values  # Subsurface runoff (out of gridcell)
        o.variables['Qsm'][:] = no_data_1D  # Snowmelt (solid to liquid)
        o.variables['Qfz'][:] = no_data_1D  # Re-freezing of water in the snow (liquid to solid)
        o.variables['DelSoilMoist'][:] = no_data_1D  # Change in soil moisture (increase)
        o.variables['DelSWE'][:] = no_data_1D  # Change in snow water equivalent (increase)
        o.variables['DelIntercept'][:] = no_data_1D  # Change in interception storage (increase)
        o.variables['Qirrig'][:] = no_data_1D  # Anthropogenic water flux from irrigation (increase)
        # Surface state variables
        o.variables['SnowT'][:] = no_data_1D  # Snow surface temperature
        o.variables['VegT'][:] = no_data_1D  # Vegetation canopy temperature
        o.variables['BaresoilT'][:] = no_data_1D  # Temperature of bare soil (skin)
        o.variables['AvgSurfT'][:] = data['AvgSurfT'].values  # Average surface temperature (skin)
        o.variables['RadT'][:] = no_data_1D  # Surface radiative temperature
        o.variables['Albedo'][:] = data['Albedo'].values  # Surface albedo
        o.variables['SWE'][:] = no_data_1D  # Snow water equivalent
        o.variables['SurfStor'][:] = no_data_1D  # Surface water storage
        o.variables['SnowFrac'][:] = no_data_1D  # Snow covered fraction
        o.variables['SAlbedo'][:] = no_data_1D  # Snow albedo
        o.variables['CAlbedo'][:] = no_data_1D  # Vegetation canopy albedo
        o.variables['UAlbedo'][:] = no_data_1D  # Urban canopy albedo
        o.variables['LAI'][:] = no_data_1D  # Leaf area index
        o.variables['RoofSurfT'][:] = no_data_1D  # Roof surface temperature (skin)
        o.variables['WallSurfT'][:] = no_data_1D  # Wall surface temperature (skin)
        o.variables['RoadSurfT'][:] = no_data_1D  # Road surface temperature (skin)
        o.variables['TairSurf'][:] = no_data_1D  # Near surface air temperature (2m)
        o.variables['TairCanyon'][:] = no_data_1D  # Air temperature in street canyon (bulk)
        o.variables['TairBuilding'][:] = no_data_1D  # Air temperature in buildings (bulk)
        # Sub-surface state variables **** TWO DIMENSIONAL ****
        o.variables['SoilMoist'][:, :] = data_SoilMoist  # Average layer soil moisture
        o.variables['SoilTemp'][:, :] = data_SoilTemp  # Average layer soil temperature
        # Evaporation components
        o.variables['TVeg'][:] = no_data_1D  # Vegetation transpiration
        o.variables['ESoil'][:] = no_data_1D  # Bare soil evaporation
        o.variables['RootMoist'][:] = no_data_1D  # Root zone soil moisture
        o.variables['SoilWet'][:] = no_data_1D  # Total soil wetness
        o.variables['ACond'][:] = no_data_1D  # Aerodynamic conductance
        # Forcing data (at forcing height)
        o.variables['SWdown'][:] = data['SWdown'].values  # Downward shortwave radiation
        o.variables['LWdown'][:] = data['LWdown'].values  # Downward longwave radiation
        o.variables['Tair'][:] = data['Tair'].values  # Air temperature
        o.variables['Qair'][:] = data['Qair'].values  # Specific humidity
        o.variables['PSurf'][:] = data['PSurf'].values  # Air pressure
        o.variables['Wind'][:] = data['Wind'].values  # Wind speed