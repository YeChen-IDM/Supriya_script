# dtk analyze -id 6c3bb92e-9b78-e711-9401-f0921c16849d -a analyze_TotalCases.py
from site_Zambia import ZambiaCalibSite

site = ZambiaCalibSite(
        analyzers = {
            'NationalPrevalenceAnalyzer': 0.74,
            'ProvincialPrevalenceAnalyzer': 0.12,
            'PopulationAnalyzer': 0.12,
            'ProvincialARTAnalyzer': 0.02
        },
        force_apply = False,
        max_sims_per_scenario = -1
    )
analyzers = site.get_analyzers()
