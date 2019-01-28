import importlib
import logging
from calibtool.CalibSite import CalibSite

logger = logging.getLogger(__name__)

class ZambiaCalibSite(CalibSite):
    metadata = {}

    def __init__(self, **kwargs):
        # required kwargs first
        analyzers = kwargs.get('analyzers', None)
        if not analyzers:
            raise Exception('An analyzer dictionary must be provided to the \'analyzers\' argument.')

        # reference_data must be supplied as a kwarg and is simply stored for direct use by analyzers
        self.reference_data = kwargs.get('reference_data', None)
        if not self.reference_data:
            raise Exception('Obs/reference data object must be provided to the \'reference_data\' argument.')

        # optional kwargs
        force_apply = kwargs.get('force_apply', False)
        max_sims_per_scenario = kwargs.get('max_sims_per_scenario', -1)

        self.analyzers = []
        for analyzer, weight in analyzers.items():
            AnalyzerClass = getattr(importlib.import_module('hiv.analyzers.%s'%analyzer), analyzer)

            self.analyzers.append(
                AnalyzerClass(
                    self,
                    weight = weight,

                    force_apply = force_apply,
                    max_sims_per_scenario = max_sims_per_scenario,

                    reference_year = 2010, #matches ingest form (2010.5 in integer form)
                    reference_population = 5907170, #corresponds to total 15-49 population in Zambia
                    age_min = 15, #census data starts at 0 - does this need to match here?
                    age_max = 49, #census data ends at 49 - same question as above

                    node_map = {
                        1: "Central",
                        2: "Copperbelt",
                        3: "Eastern",
                        4: "Luapula",
                        5: "Lusaka",
                        6: "Muchinga",
                        7: "Northern",
                        8: "Northwestern",
                        9: "Southern",
                        10: "Western"
                    },

                    basedir = '.',
                    fig_format = 'png',
                    fig_dpi = 600,
                    verbose = True
                )
            )

        # Must come at the end:
        super(ZambiaCalibSite, self).__init__('SEARCH')

    def get_setup_functions(self):
        return [ ]

    def get_analyzers(self):
        return self.analyzers
