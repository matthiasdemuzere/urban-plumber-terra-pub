This repository is developed to run **TERRA v11** in line with the the [Urban-PLUMBER](https://urban-plumber.github.io/) requirements.


TERRA?
--------

TERRA is the land-surface scheme in COSMO-CLM, the NWP and regional climate model developed by DWD. Its urban parameterization TERRA_URB uses the Semi-empirical URban canopy dependencY algorithm (SURY) to condense three dimensional urban canopy information (if available) into a limited number of bulk properties, thereby strongly reducing the computational cost without loss of performance. 

Unfortunately, the TERRA stand-alone (TSA) version used in the Urban-Plumber project is not developed alongside the developments in the TERRA version embedded in COSMO-CLM. As a result:
* this version does not reflect the most recent TERRA version as used in COSMO-CLM. Note here that a new state-of-the-art TSA is being developed, yet this version does not yet contain the urban parameterization. For more information on this new TSA version, please contact [Jean-Marie Bettems](mailto:Jean-Marie.Bettems@meteoswiss.ch), (MeteoSwiss).
* the calculation of some external parameters are not available in the way they are used in coupled model experiments (that use [EXTPAR](http://www.doi.org/10.1002/2013JD021267) to extract these parameters).

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
* Schulz JP, Vogel G. Improving the processes in the land surface scheme TERRA: Bare soil evaporation and skin temperature. Atmosphere (Basel). 2020;11(5):1-17. [DOI](http://doi.org/10.3390/atmos11050513).
* Varentsov, M., Samsonov, T., & Demuzere, M. (2020). Impact of Urban Canopy Parameters on a Megacity’s Modelled Thermal Environment. Atmosphere, 11(12), 1349. [DOI](http://doi.org/10.3390/atmos11121349). 


Configuration details
------------
* TERRA_URB is a bulk scheme, that works with an urban and natural tile. When running the stand-alone version, the model has to be run twice, one time for urban `CASE = urb` and one time for `CASE = veg`). 


* Afterwards, the results are weighted according to their fractions. Note that TERRA_URB does not simulate water bodies, so this fraction is neglected, and fractions are thus rescaled. 

* For the urban fraction, TERRA uses [SURY](https://github.com/hendrikwout/sury/blob/master/sury.py), the Semi-empirical URban canopy dependencY algorithm that is used to translate 3D building parameters into bulk values. 


* When running the natural tile, the natural weights need to be rescaled to unity, and this information is then integrated in the EXTPAR namelist settings via `plcov_max/min_const`, and `plant cover fraction`.


* TERRA works with a soiltype parameter look-up table, defined in `data_soil.f90`. The available soil types are `ice, rock, sand, sandy loam, loam, clay loam, clay, peat, sea water, sea ice, urban`, chosen via the input parameter `soiltyp_const` in the `EXTPARA` input namelist section. The appriopriate soil type is automatically set by selecting the soil type class that returns the minimum difference for the clay and sand fractions combined, evaluated against `topsoil_clay_fraction` and `topsoil_sand_fractio` provided in the sitedata. E.g. for Preston, the soiltype class is `sandy loam`, with clay / sand fractions of 0.10 / 0.65 & 0.18 / 0.72 for TERRA look-up & sitedata respectively.
 
  
* When coupled to COSMO:
    * min/max LAI is extracted using [EXTPAR](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013JD021267), via which - for every pixel - the min/max LAI is obtained by weighing land-cover specific LAI values with thier corresponding GLOBCOVER land use fractions. As this system is not available for TERRA stand-alone, these values are approximated by averaging GLOBCOVER's tree, grass and bare soil LAI values, and weighing them with the fractions provided by Urban Plumber.

  * AHF is taken into account - sourced from 
[Flanner (2009)](http://www.agu.org/pubs/crossref/2009/2008GL036465.shtml) | [version with errate fixed](http://clasp-research.engin.umich.edu/faculty/flanner/content/ppr/ppr_ener.pdf) | [AHF_2005_2.5min.nc](https://www.cgd.ucar.edu/tss/ahf/data/) - and added as an additional source of sensible heat flux added to the lowest model layer. In stand-alone mode, this routine is by default not availble. From _version 0.0.2b_ of this script chain, this feature is added, in a similar way as done within COSMO. Note that, for:
    *  **_BASELINE_**: the mean annual AHF(lat,lon) value from Flanner is taken
    *  **_DETAILED_**: the mean annual AHF value from the sitedata (`anthropogenic_heat_flux_mean`) is used.
    *  **_LCZ_**: the mean annual AHF value from the corresponding LCZ class is used.
  
  * Recent changes in [bare soil evaporation](http://www.doi.org/10.3390/atmos11050513) available and used in the official COSMO model were not available in TERRA TSA. These changes have been backported to the TSA version, and are used as a default in the Urban-Plumber experiments (`itype_evsl=4`, see below).


* For the **_BASELINE_** experiment, the following sitedata info is used:
    * latitude
    * measurement_height_above_ground
    * impervious_area_fraction
    * tree_area_fraction
    * grass_area_fraction
    * bare_soil_area_fraction
    * water_area_fraction
    * topsoil_clay_fraction
    * topsoil_sand_fraction
  

* For the **_DETAILED_** experiment, the following additional sitedata info is used:
    * roughness_length_momentum
    * average_albedo_at_midday
    * roof_area_fraction (input to SURY)
    * canyon_height_width_ratio (input to SURY)
    * building_mean_height (input to SURY)
    * anthropogenic_heat_flux_mean (input to Flanner's (2009) parametric implementation).



Simulation status
------------

### Phase I
| Site | Baseline | Detailed | LCZ** | Version | Code Release |
| --- | --- | --- | --- | --- | --- | 
| AU-Preston | Submitted | Submitted | Available | 1 | v0.0.1b |
| AU-Preston | Submitted | NA | NA | 2 | v0.0.2b0 |
| AU-Preston | Submitted | Submitted | Available | 4 | v0.4.0 |

### Phase II
| Site | Baseline | Detailed | LCZ** | Version | Code Release |
| --- | --- | --- | --- | --- | --- | 
| AU-Preston | Submitted | Submitted | Available | 5 | v0.5.0 |
| AU-SurreyHills | Submitted | Submitted | Available | 5 | v0.5.0 |
| CA-Sunset | Submitted | Submitted | Available | 5 | v0.5.0 |
| FI-Kumpula | Submitted | Submitted | Available | 5 | v0.5.0 |
| FI-Torni | Submitted | Submitted | Available | 5 | v0.5.0 |
| FR-Capitole | Submitted | Submitted | Available | 5 | v0.5.0 |
| GR-HECKOR | Submitted | Submitted | Available | 5 | v0.5.0 |
| JP-Yoyogi | Submitted | Submitted | Available | 5 | v0.5.0 |
| KR-Jungnang | Submitted | Submitted | Available | 5 | v0.5.0 |
| KR-Ochang | Submitted | Submitted | Available | 5 | v0.5.0 |
| MX-Escandon | Submitted | Submitted | Available | 5 | v0.5.0 |
| NL-Amsterdam | Submitted | Submitted | Available | 5 | v0.5.0 |
| PL-Lipowa | Submitted | Submitted | Available | 5 | v0.5.0 |
| PL-Narutowicza | Submitted | Submitted | Available | 5 | v0.5.0 |
| SG-TelokKurau | Submitted | Submitted | Available | 5 | v0.5.0 |
| UK-KingsCollege | Submitted | Submitted | Available | 5 | v0.5.0 |
| UK-Swindon | Submitted | Submitted | Available | 5 | v0.5.0 |
| US-Baltimore | Submitted | Submitted | Available | 5 | v0.5.0 |
| US-Minneapolis1 | Submitted | Submitted | Available | 5 | v0.5.0 |
| US-Minneapolis2 | Submitted | Submitted | Available | 5 | v0.5.0 |
| US-WestPhoenix | Submitted | Submitted | Available | 5 | v0.5.0 |

### Phase II (fix bugs in precip forcing + hydrology units)
| Site | Baseline | Detailed | LCZ** | Met Version | Code Release |
| --- | --- | --- | --- |------------| --- | 
| AU-Preston | Submitted | Submitted | Available | v3 | v0.7.0 |
| AU-SurreyHills | Submitted | Submitted | Available | v1 | v0.7.0 |
| CA-Sunset | Submitted | Submitted | Available | v1 | v0.7.0 |
| FI-Kumpula | Submitted | Submitted | Available | v1 | v0.7.0 |
| FI-Torni | Submitted | Submitted | Available | v1 | v0.7.0 |
| FR-Capitole | Submitted | Submitted | Available | v1 | v0.7.0 |
| GR-HECKOR | Submitted | Submitted | Available | v1 | v0.7.0 |
| JP-Yoyogi | Submitted | Submitted | Available | v1 | v0.7.0 |
| KR-Jungnang | Submitted | Submitted | Available | v1 | v0.7.0 |
| KR-Ochang | Submitted | Submitted | Available | v1 | v0.7.0 |
| MX-Escandon | Submitted | Submitted | Available | v1 | v0.7.0 |
| NL-Amsterdam | Submitted | Submitted | Available | v1 | v0.7.0 |
| PL-Lipowa | Submitted | Submitted | Available | v1 | v0.7.0 |
| PL-Narutowicza | Submitted | Submitted | Available | v1 | v0.7.0 |
| SG-TelokKurau | Submitted | Submitted | Available | v1 | v0.7.0 |
| UK-KingsCollege | Submitted | Submitted | Available | v1 | v0.7.0 |
| UK-Swindon | Submitted | Submitted | Available | v1 | v0.7.0 |
| US-Baltimore | Submitted | Submitted | Available | v1 | v0.7.0 |
| US-Minneapolis1 | Submitted | Submitted | Available | v1 | v0.7.0 |
| US-Minneapolis2 | Submitted | Submitted | Available | v1 | v0.7.0 |
| US-WestPhoenix | Submitted | Submitted | Available | v1 | v0.7.0 |


** For the online TERRA version, the [WUDAPT-to-COSMO](https://github.com/matthiasdemuzere/WUDAPT-to-COSMO) tool is available, that provides spatially explicit urban canopy parameters using the concept of [Local Climate Zones](www.wupdat.org) and their corresponding parameters values from [Stewart and Oke (2012)](http://10.1175/BAMS-D-11-00019.1) and [Stewart et al. (2014)](http://10.1002/joc.3746) (see also `LCZ_UCP_default.csv` under `tables/`). As this work is still under development ([Varentsov et al., 2020](https://www.mdpi.com/2073-4433/11/12/1349)), it is not considered as the default approach for TERRA. Yet this **_LCZ_** approach is performed alongside the **_BASELINE_** and **_DETAILED_** experiment, to further test it across a range of flux tower sites provided by this experiment. <br>
** For this LCZ approach, surface fractions are used from the **_BASELINE_** information (lcz tables have no info on grass/tree/bare soil), and all other info needed to produce bulk parameters with [SURY](https://github.com/hendrikwout/sury/blob/master/sury.py) from the LCZ look-up tables (building height, heigh-to-width ratio, roof fraction, facet (roof, wall, road) albedo / emissivity / heat capacity / heat conductivity, and the anthropogenic heat flux).



Notes
------------ 
* When forcing TERRA with meteo, the meteo input data needs to start at 00h00 and end at 23h or 23h30 (depending on the forcing interval). For now, missing times are filled with 0's (NaNs not allowed), and will be clipped from the analysis. No impact is expected due to 1) length spin-up (in case 0's occur at the beginning) or clipping at end for analysis (e.g. Preston case where time series ends at 13:00).
* The TERRA code is not open-source. One needs to be a member of the [COSMO-CLM community](https://wiki.coast.hzg.de/clmcom/) before getting access. Hence it is not included in this repository.


Changelog - Versions
------------
* **v0.7.0**: 
  * TERRA is now forced with exact rain and snow rates. This in contrast to previous versions where 1) precip = rain + snow and 2) precip was converted back to rain and snow within TERRA using the a temperature threshold of 273.15K. This was actually done twice (during input readig and temporal interpolation) and was causing inconsistencies between the forcing and what was actually used in TERRA. 

* **v0.6.0**: 
  * Fix bug in the units of Qs and Qsb (runoff) provided in the final ouput (in `create_final_df_from_grb()`). 
  * Within the TERRA src code, precip (rainf + snowf) is split twice into rain and snow forcing: one time in _read_statbin_, and again during _the temporal interpolation_ (both part of `terra_io _30 or _60 .f90`. This is a bug, and as such, this procedure is removed from _read_statbin_, and only kept during `the temporal interpolation`.

  
* **v0.5.0** (Phase II - first version):
  * Final production runs use 2 years of spin-up. That has found to be sufficient to stabilize initial soil moisture and surface temperatures.
  * Final output file name constructed as: `{TERRA_VERSION}_{SITE}_{EXPERIMENT}_{VERSION}.nc`. Eg. `TERRA_4.11_FR-Capitole_d_v5.nc`
  * Fix time setting for sites with `time_coverage_end` ending with 00:00:00 (was not taken into account).
  * Introduce a scaler to scale (lower) initial soil moisture. Used for `GR-HECKOR (2)`, `PL-Lipowa (2)`, `PL-Narutowicza (2)`,  `US-Minneapolis1 (3)`,  `US-Minneapolis2 (3)`, and `US-WestPhoenix (5)`.
  * Backport `itype_evsl = 4` from COSMO 5.0 code to TERRA TSA, allowing the use of the resistance version of bare soil evaporation. 
  * **BUG FIX** in `WIMPC_C`. Value was 1e6 times too large, leading to 0 Qle from urban surfaces.
  * **BUG FIX** in `_get_soil_type_class`. Now only sampling from true soil types: `['sand', 'sandy loam', 'loam', 'clay loam', 'clay', 'peat']`. See [here](https://github.com/matthiasdemuzere/urban-plumber-terra/commit/161ae970d4a6850d2ea14a5b45a7317e806a2ae9).
  * Full revision / restructuring of the code, making the simulation process more easy
  * Now only read info from .nc file (including sitedata). The .csv input is no longer used. 
  * LCZ class and netcdf meteo forcing version now read from `tables/site_lookup.csv`
  * The choice in simulation length is now simplified, by adding the `START_RUN` argument. More info below under **Execute the toolchain**.
  * Automatically assign soiltype class using provided clay and sand fractions.
  * Add check to make sure model compiled without error. If error, script aborts ...
  * A number of parameter flags / values were changed, to the following:
      * lvegadapt=.true.
      * lrootadapt=.false.
      * itype_root=2
      * itype_evsl=4
      * ROOTDP=2.5
  * Add Mikhail Varentsov as secundary contact in attributes of netcdf file.
  
* **v0.4.0** (Phase I): 
  * Allow half hourly I/O.
  * Add additional output variables: `SWup`, `LWup`, `Qs`, `Qsb`, `Albedo`, `AvgSurfT`, `Rainf`, `Snowf`. 
  * Fix bug in precipitation forcing units. 
  * Fix bug in AHF implentation (need to use local times instead of UTC)
  * Introduce more flexible debug system: set `debug=True` will run case for first three months.
  * Small edits to global attributes in netcdf output.
  * `dt` refers to internal model timestep. Keep this constant (60s). The final output frequency scales with `dt` and `nout_interval`: `output frequency = dt * nout_interval`, so 1800s = 60 (`dt`) * 30 (`nout_interval`).
  * Final output file name now also reflects temporal output frequency: `{TERRA_VERSION}_{SITE}_` **INPUT_INTERVAL**` _{CONFIG}{VERSION}.nc`
  

* **v0.0.2b0**: Updates based on feedback, including:
  * Using H, htw and roof_f as input to SURY for the **_DETAILED_** experiment. This affects bulk albedo, emissivity, heat capacity and conductivity and roughness length.
  * Integration of AHF (Qanth) as used in COSMO-CLM, for **_BASELINE_** (using Flanner's value), for **_DETAILED_** (using sitedata info) and for **_LCZ_** (using the LCZ class value). For all, Flanner's parameterisation is applied to prescribe a diurnal and seasonal cycle. 
 

* **v0.0.1b**: initial commit
