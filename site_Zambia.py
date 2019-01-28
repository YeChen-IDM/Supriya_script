import importlib
import logging
from calibtool.CalibSite import CalibSite
from hiv.analysis.HIVCalibSite import HIVCalibSite
logger = logging.getLogger(__name__)

class ZambiaCalibSite(CalibSite):
    metadata = {}

    def __init__(self, **kwargs):
        # required kwargs first
        analyzers = kwargs.get('analyzers', None)
        if not analyzers:
            raise Exception('An analyzer dictionary must be provided to the \'analyzers\' argument.')
        site_info = kwargs.get('site_info', None)

        # reference_data must be supplied as a kwarg and is simply stored for direct use by analyzers
        self.reference_data = kwargs.get('reference_data', None)
        if not self.reference_data:
            raise Exception('Obs/reference data object must be provided to the \'reference_data\' argument.')

        # optional kwargs
        force_apply = kwargs.get('force_apply', False)
        max_sims_per_scenario = kwargs.get('max_sims_per_scenario', -1)

        site = HIVCalibSite(analyzers=analyzers,
                            reference_data=self.reference_data,
                            site_data=site_info,
                            force_apply=force_apply,
                            max_sims_per_scenario=max_sims_per_scenario)
        self.analyzers = site.analyzers

        # Must come at the end:
        super(ZambiaCalibSite, self).__init__('SEARCH')

    def get_setup_functions(self):
        return [ ]

    def get_analyzers(self):
        return self.analyzers
