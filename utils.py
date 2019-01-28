import json, os
from dtk.utils.builders.TaggedTemplate import CampaignTemplate


def make_campaign_template(
        plugin_files_dir = os.path.join('..', 'InputFiles'),
        is_fasttrack = False,
        with_placebo = False,
        with_commercial_sex_events = True,
        with_medical_male_circumcision_events = True,
        with_traditional_male_circumcision_events = True,
        with_sti_coinfection_events = True,
        with_risk_reduction_events = True,
        with_PrEP = True,
        with_PrEP_cascade = True,
        additional_overlay_filepaths = []
    ):

    base = os.path.join(plugin_files_dir, 'campaign_Zambia_DMPA_control_group.json')
    # base = os.path.join(plugin_files_dir, 'campaign_SEARCH_Community_SQ_Stop_After_2016_2D_after_adding_baseline_linkage.json')
    overlay_event_filepaths_list = []

    return CampaignTemplate.from_files( base, overlay_event_filepaths_list)

def loadOptimConfigFromFile(opt_file):
    with open(opt_file, 'r') as json_file:
        J = json.loads(json_file.read())

    jp = []
    x = []
    for comp in J['Model']['Components']:
        print('-'*80)
        for i, (c, value) in enumerate(zip(comp['json_paths'], comp['x_current_postmap'])):
            print('%d: %s = %s' % (i, c, value))
            path = str(c);
            if '.' not in path:
                jp.append(path)
                x.append(value)
                next

            tok = path.split('.')
            if tok[0] in ['CONFIG', 'CAMPAIGN', 'DEMOGRAPHICS']:
                path = '.'.join(tok[1:])
            else:
                print("WARNING, confused about parameter: ", path)
            jp.append(path)
            x.append(value)


    return (jp,x)
