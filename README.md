
This repository is developed to run **TERRA v11** in line with the the [Urban-PLUMBER](https://urban-plumber.github.io/) requirements.  


TERRA?
--------

TERRA in the land-surface scheme in COSMO-CLM, the NWP and regional climate model developed by DWD. Its urban parameterization TERRA_URB uses the Semi-empirical URban canopy dependencY algorithm (SURY) to condense the three dimensional urban canopy information to a limited number of bulk properties, thereby strongly reducing the computational cost without loss of performance. 

Unfortunately, the TERRA stand-alone version used in this experiment is not developed alongside the developments in the TERRA version embedded in COSMO-CLM. As a result:
* some of the recent developments are not implemented (e.g. [bare soil evaporation](http://www.doi.org/10.3390/atmos11050513))
* the calculation of some external parameters are not available as they are used in coupled model experiments, as this is done with [EXTPAR](http://www.doi.org/10.1002/2013JD021267).

See below for more details on the configuration.  

**References (selected)**:

* Wouters, H.; Demuzere, M.; Ridder, K.D.; van Lipzig, N.P.M. (2015) The impact of impervious water-storage parametrization on urban climate modelling.
Urban Climate, 11, 24–50. [DOI](http://doi.org/10.1016/j.uclim.2014.11.005)
* Trusilova, K., Schubert, S., Wouters, H., Früh, B., Grossman-Clarke, S., Demuzere, M., & Becker, P. (2016). The urban land use in the COSMO-CLM model: a comparison of three parameterizations for Berlin. Meteorologische Zeitschrift, 25(2), 231–244. [DOI](https://doi.org/10.1127/metz/2015/0587)
* Wouters, H.; Demuzere, M.; Blahak, U.; Fortuniak, K.; Maiheu, B.; Camps, J.; Tielemans, D.; van Lipzig,
N.P. (2016) Efficient urban canopy parametrization for atmospheric modelling: description and application with the COSMO-CLM model (version 5.0_clm6) for a Belgian Summer. Geoscientific Model Development,
9, 3027–3054. [DOI](http://doi.org/10.5194/gmd-9-3027-2016)
* Wouters, H.; De Ridder, K.; Poelmans, L.; Willems, P.; Brouwers, J.; Hosseinzadehtalaei, P.; Tabari, H.;
Vanden Broucke, S.; van Lipzig, N.P.M.; Demuzere, M. (2017) Heat stress increase under climate change twice as large in cities as in rural areas: A study for a densely populated midlatitude maritime region. Geophysical Research Letters, 44, 8997–9007. [DOI](http://doi.org/10.1002/2017GL074889).
* Demuzere, M.; Harshan, S.; Järvi, L.; Roth, M.; Grimmond, C.S.B.; Masson, V.; Oleson, K.W.; Velasco,
E.; Wouters, H. (2017) Impact of urban canopy models and external parameters on the modelled urban energy
balance in a tropical city. Quarterly Journal of the Royal Meteorological Society, 143, 1581–1596.
[DOI](http://doi.org/10.1002/qj.3028).
* Varentsov, M., Samsonov, T., & Demuzere, M. (2020). Impact of Urban Canopy Parameters on a Megacity’s Modelled Thermal Environment. Atmosphere, 11(12), 1349. [DOI](http://doi.org/10.3390/atmos11121349).  


Instructions
------------

It is advised to use a python virtual environment:

1. Go into python environments directory: `cd ~/python-environments/`
2. Create virtual environment: `python3 -m venv up_terra` or `virtualenv --python=/usr/bin/python3.6 up_terra`
3. Go to SCRIPTDIR: `cd SCRIPTDIR`
4. Install module requirements: `~/python-environments/up_terra/bin/pip3 install -r requirements.txt`
5. For in terminal: Activate environment: `. ~/python-environments/up_terra/bin/activate`

One the environment is created and activated, execute the script chain with:
```
python3 run_chain.py ARG1 ARG2 ARG3 ARG4 ARG5

with:
  ARG1: site name (eg. AU-Preston)
  ARG2: config    ('b' or 'd' or 'lcz')
  ARG3: version   (integer)
  ARG4: input/output interval (integer, 30 or 60 [minutes])
  ARG5: debug mode, boolean (True or False). When True, only run for 1 year.
```


**Note**: The `requirements.txt` can be generated using `pipreqs`: 
```
cd /SCRIPT/DIR/
pipreqs --ignore=terra/ .
```


Simulation status
------------
| Site | Baseline | Detailed | LCZ** | Version |
| --- | --- | --- | --- | --- |
| AU-Preston | Submitted | Submitted | Available | 1 |
| AU-Preston | progress | progress | progress | 2 |

** For the online TERRA version, we are currently developing the [WUDAPT-to-COSMO](https://github.com/matthiasdemuzere/WUDAPT-to-COSMO) generator, that provides spatially explicit urban canopy parameters using the concept of [Local Climate Zones](www.wupdat.org) and their corresponding parameters values from [Stewart and Oke (2012)](http://10.1175/BAMS-D-11-00019.1) and [Stewart et al. (2014)](http://10.1002/joc.3746) (see also `LCZ_UCP_default.csv` under `tables/`). As this work is still under development ([Varentsov et al., 2020](https://www.mdpi.com/2073-4433/11/12/1349)), it is not considered as the default approach for TERRA. Yet this **_LCZ_** approach is performed alongside the **_BASELINE_** and **_DETAILED_** experiment, to further test it across a range of flux tower sites provided by this experiment. <br>
** For this LCZ approach, surface fractions are used from the **_BASELINE_** information (lcz tables have no info on grass/tree/bare soil), and all other info needed to produce bulk parameters with [SURY](https://github.com/hendrikwout/sury/blob/master/sury.py) from the LCZ look-up tables (building height, heigh-to-width ratio, roof fraction, facet (roof, wall, road) albedo / emissivity / heat capacity / heat conductivity, and the anthropogenic heat flux).   


Configuration details
------------
* TERRA_URB is a bulk scheme, that works with an urban and natural tile. 
When running the stand-alone version, the model has to be run twice, one time for urban `CASE = urb` and 
one time for `CASE = veg`). 

* Afterwards, the results are weighted according to their fractions. Note that TERRA_URB does not simulate water bodies, so this fraction is neglected, and fractions are thus rescaled. 

* When running the natural tile, the natural weights need to be rescaled to unity,
and this information is then integrated in `EXTPARA - plcov_max_const & plcov_min_const max / min plant cover fraction)`.

* When coupled to COSMO, min/max LAI is extracted using [EXTPAR](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013JD021267), in which for every pixel min/max LAI is obtained by weighing land-cover specific LAI values with thier corresponding GLOBCOVER land use fractions. As this system is not available for TERRA stand-alone, these values are approximated by averaging GLOBCOVER's tree, grass and bare soil LAI values, and weighing them with the fractions provided by Urban Plumber.

* When coupled to COSMO, AHF is taken into account - sourced from 
[Flanner (2009)](http://www.agu.org/pubs/crossref/2009/2008GL036465.shtml) | [version with errate fixed](http://clasp-research.engin.umich.edu/faculty/flanner/content/ppr/ppr_ener.pdf) | [AHF_2005_2.5min.nc](https://www.cgd.ucar.edu/tss/ahf/data/) - and added as an additional source of sensible heat flux added to the lowest model layer. In stand alone mode, this routine is by default not availble. From _version 0.0.2b_ of this script chain, this feature is added, in a similar way as done within COSMO. Note that, for:
  *  **_BASELINE_**: the mean annual AHF(lat,lon) value from Flanner is taken
  *  **_DETAILED_**: the mean annual AHF value from the sitedata (`anthropogenic_heat_flux_mean`) is used.

* For the **_BASELINE_** experiment, the following sitedata info is used:
    * latitude
    * measurement_height_above_ground
    * impervious_area_fraction
    * tree_area_fraction
    * grass_area_fraction
    * bare_soil_area_fraction
    * water_area_fraction
    
* For the **_DETAILED_** experiment, the following sitedata info is used:
    * roughness_length_momentum
    * average_albedo_at_midday
    * roof_area_fraction (input to SURY)
    * canyon_height_width_ratio (input to SURY)
    * building_mean_height (input to SURY)
    * anthropogenic_heat_flux_mean (input to Flanner's (2009) parametric implementation.
    
* Note that [SURY](https://github.com/hendrikwout/sury/blob/master/sury.py) is the Semi-empirical URban canopy dependencY algorithm, that is used to translate 3D building parameters into bulk values. 

* When forcing TERRA with meteo, the `NRDAYS` needs to be provided. As such, info for all (half) hours of a day is required, NANs are not excepted. For now, missing times are filled with 0's, and will be clipped from the analysis. No impact expected due to 1) length spin-up (in case 0's occur at the beginning) or clipping at end for analysis (e.g. Preston case where time series ends at 13:00).


Notes
------------ 
* The TERRA code is not open-source. One needs to be a member of the [COSMO-CLM community](https://wiki.coast.hzg.de/clmcom/) before getting access. Hence it is not included in this repository.

Versions
------------
* **UP NEXT**: 
  * Allow half hourly input and output, compared to only hourly in version **v0.0.1b**.
  * Small edit to global attributes in netcdf output.
  * `dt` refers to internal model timestep. Keep this constant (60s). The final output frequency scales with `dt` and `nout_interval`: `output frequency = dt * nout_interval`, so 1800s = 60 (`dt`) * 30 (`nout_interval`).
  * Final output file name now also reflects temporal output frequency: `{TERRA_VERSION}_{SITE}_` **INPUT_INTERVAL**` _{CONFIG}{VERSION}.nc` 
  
* **v0.0.2b0**: Updates based on feedback, including:
  * Using H, htw and roof_f as input to SURY for the **_DETAILED_** experiment. This affects bulk albedo, emissivity, heat capacity and conductivity and roughness length.
  * Integration of AHF, for **_BASELINE_** (using Flanner's value), for **_DETAILED_** (using sitedata info) and for **_LCZ_** (using the LCZ class value). For all, Flanner's parameterisation is applied to prescribe a diurnal and seasonal cycle. 
 

* **v0.0.1b**: initial commit




    
    