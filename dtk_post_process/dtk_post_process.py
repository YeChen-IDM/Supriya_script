#!/usr/bin/python

from __future__ import print_function

import os
import pandas as pd
import time
import math

# Relative to DTK output folder:
ByAgeAndGender_filename     = "ReportHIVByAgeAndGender.csv"
Calibration_filename        = "Calibration.json"

AGGREGATED_NODE = 0 # reserved node number for aggregated aka 'National' processing

OUTPUT_DIRECTORY = 'output'

def timing(f, message):
    if message is not None:
        print('\r' + message, end='')
    t_start = time.clock()
    result = f()
    t_end = time.clock()
    if message is not None:
        print('{0:>10f}'.format(t_end-t_start))

    return result


def application(output_dir):
    print("Hello from Python!")
    print("Started Python post processing  @ " + time.asctime())
    print("Current working directory is: " + os.getcwd())

    filename = os.path.join(output_dir, ByAgeAndGender_filename)

    if not os.path.isfile(filename):
        print("!!!! Can't open " + filename + "!")
        return

    print('Loading file: %s' % filename)
    data = timing(lambda: pd.read_csv(filename), message='Load data:    ')

    reports = get_reports(data)

    data.columns = map(lambda s: s.strip().replace(' ', '_').replace('(', '').replace(')', ''), data.columns)
    node_ids = sorted([int(node_id) for node_id in data.NodeId.unique()])

    post_process_dir = 'post_process'
    directory = os.path.join(output_dir, post_process_dir)
    if not os.path.exists(directory):
        os.makedirs(directory)

    for report in reports:
        result = timing(lambda: process_report(report, data, node_ids), message=report['Name'])
        result.to_csv(os.path.join(output_dir, post_process_dir, '%s.csv' % report['Name']))

    print("Finished Python post processing @ " + time.asctime())


def get_reports(data):
    entries = []

    first_year = int(math.ceil(data.Year.min()))
    last_prevalence_year = int(data.Year.max())

    entry = {'Name': 'Prevalence',
             'Type': 'Prevalence',
             'Year': [2002.5, 2007.5, 2014.5, 2018.5],
             'AgeBins': [(15, 50)],
             'Gender': [0, 1],
             'ByNode': 'Both',
             'Map': lambda rows: rows.sum(),
             'Reduce': lambda row: float(row.Infected) / float(row.Population) if row.Population > 0 else 0}
    entries.append(entry)

    entry = {'Name': 'Population',
             'Type': 'Prevalence',
             'Year': [2010.5, 2015.5],
             'AgeBins': [(0, 5), (5, 10), (10, 15), (15, 20), (20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 50), (50, 55), (55, 60), (60, 65), (65, 70), (70, 75), (75, 80), (80, 100), (80, 85), (85, 90), (90, 100), (0, 100)],
             'Gender': [0, 1],
             'ByNode': True,
             'Map': lambda rows: rows.sum(),
             'Reduce': lambda row: row.Population}
    entries.append(entry)

    entry = {'Name': 'OnART',
             'Type': 'Prevalence',
             'Year': [2014.5, 2015.5, 2016.5, 2017.5],
             'AgeBins': [(0,15),(15,100)],
             'Gender': [0, 1],
             'ByNode': True,
             'Map': lambda rows: rows.sum(),
             'Reduce': lambda row: row.On_ART}
    entries.append(entry)

    # entry = {'Name': 'Incidence',
    #          'Type': 'Incidence',
    #          'Year': list(range(first_year, last_prevalence_year + 1)),  # must be integers
    #          'AgeBins': [(15, 25), (15, 50)],
    #          'Gender': [0, 1],
    #          'ByNode': True,
    #          'Map': lambda rows: rows.sum(),
    #          'Reduce': compute_incidence}
    # entries.append(entry)

    return entries


def compute_incidence(row):
    newly_infected_annualized = float(row.newly_infected_annualized)
    population = float(row.Population)
    newly_infected_mid = float(row['Newly_Infected'])
    if population > 0:
        national_incidence = newly_infected_annualized / (population - row.Infected - newly_infected_mid)
    else:
        national_incidence = 0
    # print('%s / (%s - %s) = %s' % (newly_infected_annualized, population, newly_infected_mid, national_incidence))
    # print(row)
    return national_incidence


def add_year_in(df):
    """
    Determine the year a particular date/timeframe occurred in (assumes tail-end dates).
    For use in 'E'-type sheets.
    """
    import numpy as np
    year_in = np.ceil(df['Year']) - 1
    return df.assign(year_in=year_in)


def preprocess_for_incidence(all_data):
    input_stratifiers = ['Year', 'NodeId', 'Gender', 'Age', 'IsCircumcised']
    grouping_stratifiers = ['NodeId', 'Gender', 'Age', 'year_in', 'IsCircumcised']

    # add the year each row is in
    data = add_year_in(all_data)
    data = data[list(set(input_stratifiers + ['year_in', 'Year', 'Newly_Infected', 'Infected', 'Population']))]

    # yearly incidence count, reported on 0.5 year to line up with the data it will be combined with in a calculation
    summed = data.drop('Year', axis=1)
    summed = summed.groupby(grouping_stratifiers).sum().reset_index()[grouping_stratifiers + ['Newly_Infected']]
    summed = summed.assign(year_in=summed['year_in'] + 0.5)
    summed = summed.rename(columns={'year_in': 'Year', 'Newly_Infected': 'newly_infected_annualized'})

    # merge into original dataframe
    summed = summed.set_index(input_stratifiers)
    data = data.set_index(input_stratifiers)
    data = data.merge(summed, how='inner', left_index=True, right_index=True).reset_index()

    # convert back to integer year; this is needed to report results on integer years as requested
    data = data.drop('Year', axis=1)
    data = data.rename(columns={'year_in': 'Year'})

    # drop circum. column; only because we added it
    data = data.drop('IsCircumcised', axis=1)
    return data


def get_blank_dataframe():
    output_stratifiers = ['Year', 'Node', 'Gender', 'AgeBin']
    output = pd.DataFrame(columns=(output_stratifiers + ['Result']))
    return output


def process_nodes(data, output, year, gender, min_age, max_age, report, node_ids):
    gender_str = 'Male' if gender == 0 else 'Female'
    for node_id in node_ids:
        rows = timing(lambda: data[
            (data.Year == year) &
            (data.Gender == gender) &
            (data.NodeId == node_id) &
            (data.Age >= min_age) &
            (data.Age < max_age)], message=None)  # 'Select rows: ')
        mapping = timing(lambda: report['Map'](rows), message=None)  # 'Sum data: ')
        try:
            result = report['Reduce'](mapping)
        except AttributeError:
            print(' -- FAILED!')
            return output
        output.loc[output.shape[0], :] = (year, node_id, gender_str, '[%d:%d)' % (min_age, max_age), result)
    return None


def process_aggregated_nodes(data, output, year, gender, min_age, max_age, report):
    gender_str = 'Male' if gender == 0 else 'Female'
    rows = timing(
        lambda: data[(data.Year == year) & (data.Gender == gender) & (data.Age >= min_age) & (data.Age < max_age)],
        message=None)  # 'Select rows: ')
    mapping = timing(lambda: report['Map'](rows), message=None)
    try:
        result = report['Reduce'](mapping)
    except AttributeError:
        print(' -- FAILED!')
        return output
    output.loc[output.shape[0], :] = (year, AGGREGATED_NODE, gender_str, '[%d:%d)' % (min_age, max_age), result)
    return None


def process_report(report, all_data, node_ids):
    # preprocessing for incidence types
    if report['Type'] == 'Incidence':
        data = preprocess_for_incidence(all_data)
    else:
        data = all_data

    # determine which type(s) of processing to do
    process_by_node = True if report['ByNode'] in [True, 'Both'] else False
    process_aggregated = True if report['ByNode'] in [False, 'Both'] else False

    # Create output data
    output = get_blank_dataframe()
    for year in report['Year']:
        for gender in report['Gender']:
            for min_age, max_age in report['AgeBins']:
                if process_by_node:
                    process_nodes(data, output, year, gender, min_age, max_age, report, node_ids)

                if process_aggregated:
                    process_aggregated_nodes(data, output, year, gender, min_age, max_age, report)

    output.set_index(['Year', 'Node', 'Gender', 'AgeBin'], inplace=True)
    return output


if __name__ == '__main__':
    application(OUTPUT_DIRECTORY)
