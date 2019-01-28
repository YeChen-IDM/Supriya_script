import copy
import math
import numpy as np
import os
import random
from scipy.special import gammaln # for calculation of mu_r

from calibtool.algorithms.OptimTool import OptimTool
from calibtool.CalibManager import CalibManager
from calibtool.plotters.LikelihoodPlotter import LikelihoodPlotter
from calibtool.plotters.OptimToolPlotter import OptimToolPlotter
from calibtool.resamplers.CramerRaoResampler import CramerRaoResampler
from calibtool.resamplers.RandomPerturbationResampler import RandomPerturbationResampler
from dtk.utils.builders.ConfigTemplate import ConfigTemplate
from dtk.utils.builders.TemplateHelper import TemplateHelper
from dtk.utils.builders.TaggedTemplate import DemographicsTemplate
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.utils.observations.utils import parse_ingest_data_from_xlsm
from simtools.SetupParser import SetupParser
from hiv.analysis.HIVCalibSite import HIVCalibSite

# local directory files
from utils import make_campaign_template


SetupParser.default_block= "HPC"

BASE_POPULATION_SCALE_FACTOR = 0.002
N_REPLICATES = 3

#
#  ***************************** PRINCIPAL SECTION OF THIS EXAMPLE ********************************
#

# The excel file with parameter, analyzer, and reference data to parse
ingest_xlsm_filename = os.path.join('Data', 'Zambia_calibration_ingest_form_v3.xlsm')

# params is a dict like below, reference is a PopulationObs object, analyzers is a dictionary of AnalyzerClass: WEIGHT
params, site_info, reference, analyzers = parse_ingest_data_from_xlsm(filename=ingest_xlsm_filename)

reference.fix_age_bins() # currently necessary, converts from [0:5) -> [0, 5) AgeBin format

# ck4, buggy, but may be necessary re: prior calibration
# Restrict to only those params that are initially active
# params = [p for p in params if p['Dynamic']]

sites = [HIVCalibSite(
            analyzers=analyzers,
            site_data=site_info,
            reference_data=reference,
            force_apply = True # Set to True to re-download and call "apply" again
            )
        ]

#
# ***************************** END PRINCIPAL SECTION OF THIS EXAMPLE ****************************
#

# dtk analyze compatibility
site = sites[0]
analyzers = site.analyzers

reserved_characters = ['<', '>', ':', '"', '/', '\\', '|', '?', '*']
for analyzer in analyzers:
    for c in reserved_characters:
        if c in analyzer.uid:
            analyzer.uid = analyzer.uid.replace(c, '-')

plotters = [
    LikelihoodPlotter(),
    OptimToolPlotter()
]

# Setting up our model configuration from templates

dir_path = os.path.dirname(os.path.realpath(__file__))
template_files_dir = os.path.join(dir_path, 'InputFiles', 'Templates')
static_files_dir = os.path.join(dir_path, 'InputFiles', 'Static')

demog = DemographicsTemplate.from_file( os.path.join(static_files_dir,'Demographics.json') )
demog_pfa = DemographicsTemplate.from_file( os.path.join(template_files_dir, 'PFA_Overlay.json') )
demog_acc = DemographicsTemplate.from_file( os.path.join(template_files_dir, 'Accessibility_and_Risk_IP_Overlay.json') )
demog_asrt = DemographicsTemplate.from_file( os.path.join(template_files_dir, 'Risk_Assortivity_Overlay.json') )
cfg = ConfigTemplate.from_file(os.path.join(template_files_dir, 'config.json'))
cpn = make_campaign_template(template_files_dir)

# For quick test simulations, this is set to a very low value
static_params = {'Base_Population_Scale_Factor': BASE_POPULATION_SCALE_FACTOR} #0.025}
cfg.set_params(static_params)

# Prepare templates
templates = TemplateHelper()
table_base = {
    'ACTIVE_TEMPLATES': [cfg, cpn, demog_pfa, demog_acc, demog_asrt],
    'TAGS': {'Scenario':'StatusQuo_Baseline', 'pyOptimTool':None}
}

config_builder = DTKConfigBuilder()
config_builder.ignore_missing = True

def constrain_sample( sample ):
    if 'Pr Ex Trns Male LOW' and 'Pr Ex Trns Male MED' in sample:
        sample['Pr Ex Trns Male LOW'] = min( [sample['Pr Ex Trns Male LOW'], sample['Pr Ex Trns Male MED']] )
    if 'Pr Ex Trns Fem LOW' and 'Pr Ex Trns Fem MED' in sample:
        sample['Pr Ex Trns Fem LOW'] = min( [sample['Pr Ex Trns Fem LOW'], sample['Pr Ex Trns Fem MED']] )

    return sample

def map_sample_to_model_input(config_builder, s):
    table = copy.deepcopy(table_base)
    table['Run_Number'] = random.randint(0, 65535) # Random random number seed

    sample = copy.deepcopy(s)

    if 'LOG Base Infectivity' in sample:
        value = sample.pop('LOG Base Infectivity')
        table['Base_Infectivity'] = np.exp(value)

    if 'PreART Link Min' and 'PreART Link Max' in sample:
        min_value = sample.pop('PreART Link Min')
        max_value = sample.pop('PreART Link Max')
        if max_value > min_value:
            table['Actual_IndividualIntervention_Config__KP_PreART_Link.Ramp_Min'] = min_value
            table['Actual_IndividualIntervention_Config__KP_PreART_Link.Ramp_Max'] = max_value
        else:
            table['Actual_IndividualIntervention_Config__KP_PreART_Link.Ramp_Min'] = max_value
            table['Actual_IndividualIntervention_Config__KP_PreART_Link.Ramp_Max'] = min_value

    if 'Male To Female Young' and 'Male To Female Old' in sample:
        young = sample.pop('Male To Female Young')
        old = sample.pop('Male To Female Old')
        table['Male_To_Female_Relative_Infectivity_Multipliers'] = [young, young, old]

    risk_reduction_fraction = sample.pop( 'Risk Reduction Fraction' ) if 'Risk Reduction Fraction' in sample else float('NaN')
    risk_ramp_rate = sample.pop( 'Risk Ramp Rate' ) if 'Risk Ramp Rate' in sample else float('NaN')
    risk_ramp_midyear = sample.pop( 'Risk Ramp MidYear' ) if 'Risk Ramp MidYear' in sample else float('NaN')

    for province in ['Central', 'Copperbelt', 'Eastern', 'Luapula', 'Lusaka', 'Muchinga', 'Northern', 'Northwestern', 'Southern', 'Western']:
        key = '%s: LOW Risk' % province
        if key in sample:
            value = sample.pop(key)
            param = 'Initial_Distribution__KP_Risk_%s' % province
            table[param] = [value, 1-value, 0]

        if not math.isnan(risk_reduction_fraction):
            param = 'Actual_IndividualIntervention_Config__KP_Medium_Risk_%s.Ramp_Max' % province
            table[param] = risk_reduction_fraction

        if not math.isnan(risk_ramp_rate):
            param = 'Actual_IndividualIntervention_Config__KP_Medium_Risk_%s.Ramp_Rate' % province
            table[param] = risk_ramp_rate

        if not math.isnan(risk_ramp_midyear):
            param = 'Actual_IndividualIntervention_Config__KP_Medium_Risk_%s.Ramp_MidYear' % province
            table[param] = risk_ramp_midyear

    if 'Risk Assortivity' in sample:
        v = sample.pop('Risk Assortivity')
        table['Weighting_Matrix_RowMale_ColumnFemale__KP_RiskAssortivity'] = [
            [ v, 1-v, 0 ],
            [ 1-v, v, v ],
            [ 0, v, 1-v ] ]

    for p in params:
        if 'MapTo' in p:
            # try:
            value = sample.pop( p['Name'] )
            if p['Name'] == "Run_Number":
                value = int(value) # Run_Number must be integer
            # except KeyError as e:
            #     print('\n************************** couldnt find key in sample: %s\n' % str(sample))
            #     raise e
            if isinstance( p['MapTo'], list):
                for mapto in p['MapTo']:
                    table[mapto] = value
            else:
                table[p['MapTo']] = value


    for name,value in sample.items():
        print('UNUSED PARAMETER:', name)

    assert (len(sample) == 0)  # All params used
    # assert( len(sample) == 1 ) # All params used except Run_Number

    return templates.mod_dynamic_parameters(config_builder, table)

# Compute hypersphere radius as a function of the number of dynamic parameters
volume_fraction = 0.1   # fraction of N-sphere area to unit cube area for numerical derivative (automatic radius scaling with N)
num_params = len( [p for p in params if p['Dynamic']] )
r = math.exp( math.log(volume_fraction) + 1/float(num_params)*(  gammaln(num_params/2.+1) - num_params/2.*math.log(math.pi) ) )
# Check, here's the formula for the volume of a N-sphere
computed_volume_fraction = math.exp( num_params/2.* math.log( math.pi ) - gammaln(num_params/2.+1) + num_params*math.log(r))

optimtool = OptimTool(params,
    constrain_sample, # <-- Will not be saved in iteration state
    #mu_r = r,      # <-- Mean percent of parameter range for numerical derivatve.  CAREFUL with integer parameters!
    mu_r = r,
    sigma_r = r/10.,   # <-- stdev of above
    samples_per_iteration = 10, # 700 is real size, 10 is testing # 3 * 140 = 420
    center_repeats = 2, # 10 is real size, 2 is testing
    rsquared_thresh = 0.81  # <-- Linear regression goodness of fit threshold.  Above this, regression is used.  Below, use best point.
)

calib_manager = CalibManager(
    name            = 'Zambia--%s--rep%s--FS' % (BASE_POPULATION_SCALE_FACTOR, N_REPLICATES), # ZimbabweOptimization_Testing ZimbabweOptimization_IncreasedSampling ZimbabweOptimization
    config_builder  = config_builder,
    map_sample_to_model_input_fn = map_sample_to_model_input,
    sites           = sites,
    next_point      = optimtool,
    sim_runs_per_param_set = N_REPLICATES, # <-- Replicates, none needed for example
    max_iterations  = 2, # limited for example
    plotters        = plotters

)

# *******************************************************************
# Resampling specific code
# *******************************************************************

# Define the resamplers to run (one or more) in list order.
resample_steps = [
    # can pass kwargs directly to the underlying resampling routines if needed
    RandomPerturbationResampler(M=60, N=2, n=1),
    CramerRaoResampler(num_of_pts=200)
]

# REQUIRED variable name: run_calib_args . Required key: 'resamplers'
run_calib_args = {
    'resamplers': resample_steps,
    'calib_manager': calib_manager
}

if __name__ == "__main__":
    SetupParser.init()
    # calib_manager.run_calibration()
    # Run the resampling
    from calibtool.ResampleManager import ResampleManager

    resample_manager = ResampleManager(steps=run_calib_args['resamplers'], calibration_manager=calib_manager)
    resample_manager.resample_and_run()
