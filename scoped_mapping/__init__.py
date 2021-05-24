from datetime import datetime
from pkg_resources import get_distribution, DistributionNotFound
from strsimpy.cosine import Cosine
import pandas as pd
import re
import requests as requests
import sqlite3
import string
import urllib
import yaml
from xml.etree import ElementTree
from tdda import rexpy

# import numpy as np
# import os
# import sys
# import rexpy

# TODO add logging
# TODO makefile for obtaining the Biosample harmonized_table.db
# TODO index sqlite tables

# ----

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    pass  # package is not installed

# ----

ols_cols = ['id', 'iri', 'short_form', 'obo_id', 'label', 'description',
            'ontology_name', 'ontology_prefix', 'type', 'is_defining_ontology']

retained_reordered_ols_cols = ['category', 'raw', 'query', 'label', 'iri', 'obo_id',
                               'ontology_name', 'ontology_prefix', 'search_rank']

ols_annotation_cols = ['iri', 'name', 'scope', 'type']

merged_cols = ['category', 'raw', 'query', 'name', 'cosine_rank', 'cosine_dist', 'obo_id', 'label',
               'search_rank', 'ontology_prefix', 'scope', 'type', 'iri', 'ontology_name']

standard_replacement_chars = '._\\- '

standard_cat_name = 'unknown'

ols_session = requests.Session()


def one_ols_submission(raw_submission, ontology_phrase='', qf_phrase='', row_req=5):
    """Get terms from OLS Search based on raw_submission string input.
    ontology_phrase = comma-delimited list of ontology names (abbreviations?) present in OLS
    qf_phrase = comma-delimited list of OLS fields to search
    row_req = the number of rows/candidate matches to request from OLS
                if there are fewer than row_req candidate matches,
                then that smaller number will be returned
                if no candidate matches are found, one blank row will be returned
    returns a dataframe
    """
    global ols_session
    userdata = {'q': raw_submission}
    if len(ontology_phrase) > 0:
        userdata['ontology'] = ontology_phrase.lower()
    if len(qf_phrase) > 0:
        userdata['queryFields'] = qf_phrase
    if row_req > 0:
        userdata['rows'] = qf_phrase
    userdata['type'] = 'class'
    userdata['exact'] = 'false'
    userdata['local'] = 'false'
    userdata['rows'] = row_req
    # userdata = {'exact': 'false',
    #              'local': 'false',
    #              'ontology': 'ncbitaxon',
    #              'q': 'homo sapiens',
    #              'rows': 5,
    #              'type': 'class',
    #             'queryFields': 'label'}
    resp = ols_session.get('https://www.ebi.ac.uk/ols/api/search', params=userdata)
    rj = resp.json()
    rjf = pd.DataFrame(rj['response']['docs'])
    rjf_rows = len(rjf.index)
    if rjf_rows < 1:
        col_count = len(ols_cols)
        blank_row = pd.Series([''] * col_count)
        rjf = rjf.append(blank_row, ignore_index=True)
        rjf.columns = ols_cols
    return rjf


def search_over_category(category_name, cat_values, ontology_phrase='', qf_phrase='', row_req=5):
    """Call one_ols_submission for multiple terms that come from some common category
        ontology_phrase = comma-delimited list of ontology names (abbreviations?) present in OLS
        qf_phrase = comma-delimited list of OLS fields to search
        row_req = the number of rows/candidate matches to request from OLS
            if there are fewer than row_req candidate matches,
            then that smaller number will be returned
            if no candidate matches are found, one blank row will be returned
        returns a list of dataframes
        """
    appended_data = []
    for one_value in cat_values:
        one_ols_result = one_ols_submission(raw_submission=one_value, ontology_phrase=ontology_phrase,
                                            qf_phrase=qf_phrase, row_req=row_req)
        one_ols_result['category'] = category_name
        one_ols_result['query'] = one_value
        one_ols_result['search_rank'] = list(range(1, len(one_ols_result.index) + 1))
        appended_data.append(one_ols_result)
    return appended_data


# had tried to escape special characters
#   that works, but raises a warning
# PEP 8: W605 invalid escape sequence '\-'
# optionally remove punctuation in separate function?
# trim whitespace?
# lowercase? OLS probably doesn't care
# important for string similarity
def whiteout(raw_string, char_string=standard_replacement_chars):
    """replaces occurrences of any character from replaced_chars
    found in raw_string with a single whitespace
    also lowercases, strips trailing or leading whitespace and collapses multiple whitespaces"""
    tidied = re.sub(r'[' + char_string + ']+', ' ', raw_string)
    tidied = re.sub(r' +', ' ', tidied)
    tidied = tidied.strip()
    tidied = tidied.lower()
    return tidied


def get_whiteout_frame(raw_list, replaced_chars=standard_replacement_chars):
    """take a list of strings and apply whiteout()
        returning a dataframe with two columns:
        the original inputs and the replaced strings"""
    wo_dict = {}
    for one_raw in raw_list:
        one_woed = whiteout(one_raw, char_string=replaced_chars)
        wo_dict[one_raw] = one_woed
    wo_frame = pd.DataFrame.from_dict(wo_dict, orient='index', columns=['woed'])
    wo_frame['raw'] = wo_frame.index
    wo_frame = wo_frame[['raw', 'woed']]
    return wo_frame


def get_wo_list(woed_frame):
    woed_list = list(set(woed_frame['woed']))
    woed_list.sort()
    return woed_list


def get_label_like(current_row):
    """Get labels, OBO synonyms (and possibly more) from submitted ontology/IRI pairs"""
    global ols_session
    anything_useful = []
    once = urllib.parse.quote(current_row['iri'], safe='')
    twice = urllib.parse.quote(once, safe='')
    term_retr_base = 'https://www.ebi.ac.uk/ols/api/ontologies/'
    term_retr_assembled = term_retr_base + current_row['ontology_name'] + '/terms/' + twice
    term_details = ols_session.get(term_retr_assembled)
    term_json = term_details.json()
    if term_json['label'] is not None:
        label_dict = {'name': term_json['label'], 'type': 'label', 'scope': 'label'}
        anything_useful.append(label_dict)
    if term_json['synonyms'] is not None:
        slistlc = [x.lower() for x in term_json['synonyms']]
        slistlc = list(set(slistlc))
        for onesyn in slistlc:
            syn_dict = {'name': onesyn, 'type': 'synonym', 'scope': 'synonym'}
            anything_useful.append(syn_dict)
    if term_json['obo_synonym'] is not None:
        anything_useful.extend(term_json['obo_synonym'])
    anno_keys = list(term_json['annotation'].keys())
    ak_lc = [one_key.lower() for one_key in anno_keys]
    syn_pattern = 'syn'
    syn_like_keys = [one_key for one_key in ak_lc if syn_pattern in one_key]
    if syn_like_keys is not None:
        for one_syn_like_key in syn_like_keys:
            vals_for_key = term_json['annotation'][one_syn_like_key]
            for one_val in vals_for_key:
                syn_dict = {'name': one_val, 'type': 'annotation', 'scope': one_syn_like_key}
                anything_useful.append(syn_dict)
    auframe = pd.DataFrame(anything_useful)
    auframe['iri'] = current_row['iri']
    return auframe


def get_csr_frame(raw_list, bad_chars=standard_replacement_chars, category_name=standard_cat_name,
                  ontoprefix='', query_fields='', rows_requested=5):
    wo_frame = get_whiteout_frame(raw_list, bad_chars)
    # get_wo_list returns a unique list
    woed_list = get_wo_list(wo_frame)
    category_search_results_list = search_over_category(category_name, woed_list, ontology_phrase=ontoprefix,
                                                        qf_phrase=query_fields, row_req=rows_requested)
    category_search_results_frame = pd.concat(category_search_results_list)
    category_search_results_frame = category_search_results_frame.merge(wo_frame,
                                                                        left_on='query', right_on='woed', how='outer')
    category_search_results_frame = category_search_results_frame[retained_reordered_ols_cols]
    return category_search_results_frame


def prep_for_label_like(category_search_results_frame):
    csr_unique = category_search_results_frame[['ontology_name', 'iri']]
    csr_unique = csr_unique[csr_unique['iri'] != '']
    csr_unique = csr_unique.drop_duplicates()
    return csr_unique


def get_bulk_label_like(label_like_prepped):
    per_iri_annotations_list = []
    for inx in label_like_prepped.index:
        one_iri_frame = get_label_like(
            {'ontology_name': label_like_prepped['ontology_name'][inx], 'iri': label_like_prepped['iri'][inx]})
        per_iri_annotations_list.append(one_iri_frame)
    all_iris_annotations_frame = pd.concat(per_iri_annotations_list)
    all_iris_annotations_frame = all_iris_annotations_frame[ols_annotation_cols]
    return all_iris_annotations_frame


# add cosine-wise ranking
# use search rank, cosine rank, type and scope for filtering/prioritizing mappings
def merge_and_compare(category_search_results_frame, bulk_label_like):
    search_annotations_merge = category_search_results_frame.merge(bulk_label_like,
                                                                   how='left', left_on='iri', right_on='iri')
    cosine_obj = Cosine(1)
    search_annotations_merge['name'] = search_annotations_merge['name'].fillna('')
    search_annotations_merge['cosine_dist'] = \
        search_annotations_merge.apply(lambda sam_row:
                                       cosine_obj.distance(sam_row.query.strip().lower(),
                                                           sam_row['name'].strip().lower()), axis=1)
    search_annotations_merge.cosine_dist = search_annotations_merge.cosine_dist.map('{:,.3f}'.format)
    all_queries = search_annotations_merge['query']
    unique_queries = list(set(all_queries))
    unique_queries.sort()
    reranked = []
    for unique_query in unique_queries:
        per_query = search_annotations_merge[search_annotations_merge['query'] == unique_query]
        pqs = per_query.sort_values(by='cosine_dist', ascending=True)
        cosine_ranks = list(range(1, len(pqs.index) + 1))
        pqs['cosine_rank'] = cosine_ranks
        reranked.append(pqs)
    reranked = pd.concat(reranked)
    return reranked


def search_get_annotations_wrapper(raw_list, bad_chars=standard_replacement_chars, cat_name=standard_cat_name,
                                   ontoprefix='', query_fields='', rr=5):
    my_category_search_results_frame = get_csr_frame(raw_list, bad_chars, category_name=cat_name,
                                                     ontoprefix=ontoprefix, query_fields=query_fields,
                                                     rows_requested=rr)
    # prep_for_label_like returns unique combinations of ontologies and term IRIs
    my_label_like_prepped = prep_for_label_like(my_category_search_results_frame)
    my_bulk_label_like = get_bulk_label_like(my_label_like_prepped)
    merged_and_compared = merge_and_compare(my_category_search_results_frame, my_bulk_label_like)
    merged_and_compared = merged_and_compared[merged_cols]
    return merged_and_compared


def read_yaml_model(modelfile):
    with open(modelfile) as file:
        inner_inferred_model = yaml.load(file, Loader=yaml.FullLoader)
    return inner_inferred_model


def get_avaialbe_enums(model):
    avaialble_enums = list(model['enums'].keys())
    avaialble_enums.sort()
    return avaialble_enums


def get_permissible_values(model, enum):
    epv = list(model['enums'][enum]['permissible_values'].keys())
    return epv


def get_best_acceptable(mappings, max_cosine=0.05):
    mappings['cosine_dist'] = mappings['cosine_dist'].astype(float)
    under_max_flag = mappings['cosine_dist'] <= max_cosine
    best_acceptable = mappings[under_max_flag]
    ba_raws = list(set(best_acceptable['raw']))
    ba_raws.sort()
    accepted_best = []
    for raw in ba_raws:
        working = best_acceptable[best_acceptable['raw'] == raw]
        min_cosine_rank = working['cosine_rank'].min()
        working = working[working['cosine_rank'] == min_cosine_rank]
        min_search_rank = working['search_rank'].min()
        working = working[working['search_rank'] == min_search_rank]
        accepted_best.append(working)
    catted_best = pd.concat(accepted_best)
    return catted_best


# # one category or enum at a time
# # has the advantage of specifying category-specific target ontologies
# # but what about the performance boost of shared queries?

def map_from_yaml(model, selected_enum, print_enums=False, bad_chars=standard_replacement_chars,
                  cat_name=standard_cat_name, ontoprefix='', query_fields=''):
    if print_enums:
        print(get_avaialbe_enums(model))
    enum_permissible_values = get_permissible_values(model, selected_enum)
    searchres_annotations = search_get_annotations_wrapper(
        enum_permissible_values,
        bad_chars=bad_chars, cat_name=cat_name, ontoprefix=ontoprefix, query_fields=query_fields)
    return searchres_annotations


# def get_no_mappings
def get_no_acceptable_mappings(all_mappings, best_acceptables):
    best_accepted_raws = set(best_acceptables['raw'])
    all_raws = set(all_mappings['raw'])
    failure_raws = all_raws - best_accepted_raws
    frl = list(failure_raws)
    failure_flag = all_mappings['raw'].isin(frl)
    failures = all_mappings[failure_flag]
    return failures


def rewrite_yaml(model, enum, best_acceptable):
    for ry_row in best_acceptable.itertuples(index=True, name='Pandas'):
        model['enums'][enum]['permissible_values'][ry_row.raw]['meaning'] = ry_row.obo_id
        model['enums'][enum]['permissible_values'][ry_row.raw]['description'] = ry_row.label


def get_sqlite_con(sqlite_file):
    sqlite_con = sqlite3.connect(sqlite_file)
    return sqlite_con


def get_package_dictionary():
    bio_s_columns = ['Name', 'DisplayName', 'ShortName', 'EnvPackage', 'EnvPackageDisplay', 'NotAppropriateFor',
                     'Description', 'Example']
    bio_s_df = pd.DataFrame(columns=bio_s_columns)
    bio_s_xml = requests.get('https://www.ncbi.nlm.nih.gov/biosample/docs/packages/?format=xml', allow_redirects=True)
    bio_s_root = ElementTree.fromstring(bio_s_xml.content)
    for node in bio_s_root:
        rowdict = {}
        for framecol in bio_s_columns:
            temp = node.find(framecol).text
            temptext = ''
            if temp is not None:
                temptext = temp
            rowdict[framecol] = temptext
        bio_s_df = bio_s_df.append(rowdict, ignore_index=True)
    return bio_s_df


def timed_query(query, connection, print_timing=False):
    start_time = datetime.now()
    if print_timing:
        print(start_time)
    result = pd.read_sql(query, connection)
    end_time = datetime.now()
    duration = end_time - start_time
    if print_timing:
        print(end_time)
        print(duration)
    return [result, duration]


def make_tidy_col(data_frame, col_in, col_out):
    punct_regex = re.compile('[%s]' % re.escape(string.punctuation))
    data_frame[col_out] = data_frame[col_in].str.lower()
    data_frame[col_out] = data_frame[col_out].str.replace(punct_regex, ' ', regex=True)
    data_frame[col_out] = data_frame[col_out].str.replace(' +', ' ', regex=True)
    data_frame[col_out] = data_frame[col_out].str.strip()
    return data_frame


def add_unique_to_list(uniquelist, non_unique):
    a = list(set(list(uniquelist)))
    b = list(set(list(non_unique)))
    c = a + b
    c = list(set(list(c)))
    c.sort()
    return c


def discover_id_pattern(example_ids):
    extracted = rexpy.extract(example_ids)
    extracted = extracted[0]
    extracted = re.sub('^\\^', '', extracted)
    extracted = re.sub('\\$$', '', extracted)
    return extracted


def add_prefix_col(dataframe, col_with_prefixes, prefix_col):
    term_list = dataframe[col_with_prefixes].str.split(':')
    prefixes = [item[0] for item in term_list]
    dataframe[prefix_col] = prefixes
    return dataframe


# def get_multi_term_patterns(dataframe, col_with_prefixes, prefix_col):
#     prefixes = dataframe[prefix_col]
#     unique_prefixes = list(set(prefixes))
#     unique_prefixes.sort()
#     inner_id_patterns = {}
#     for one_unique in unique_prefixes:
#         temp = list(ids_from_envo[col_with_prefixes][ids_from_envo[prefix_col] == one_unique])
#         extracted = rexpy.extract(temp)
#         extracted = extracted[0]
#         extracted = re.sub('^\\^', '', extracted)
#         extracted = re.sub('\\$$', '', extracted)
#         inner_id_patterns[one_unique] = extracted
#     return inner_id_patterns


def get_multi_term_patterns(dataframe, col_with_prefixes, prefix_col):
    prefixes = dataframe[prefix_col]
    unique_prefixes = list(set(prefixes))
    unique_prefixes.sort()
    inner_id_patterns = {}
    for one_unique in unique_prefixes:
        temp = list(dataframe[col_with_prefixes][dataframe[prefix_col] == one_unique])
        extracted = rexpy.extract(temp)
        extracted = extracted[0]
        extracted = re.sub('^\\^', '', extracted)
        extracted = re.sub('\\$$', '', extracted)
        inner_id_patterns[one_unique] = extracted
    return inner_id_patterns


def decompose_series(series_to_decompose, id_pattern):
    extracts = series_to_decompose.to_frame()
    extracts.columns = ['string']
    for_capture = '(' + id_pattern + ')'
    p = re.compile(for_capture)
    extracts['extract'] = extracts['string'].str.extract(p)
    extracts = extracts.fillna('')
    for_replacement = '\\[?' + id_pattern + '\\]?'
    extracts['remaining_string'] = extracts['string'].str.replace(for_replacement, '', regex=True)
    extracts = make_tidy_col(extracts, 'remaining_string', 'remaining_tidied')
    return extracts


# def extract_with_pattern(dataframe, col_to_extract, pattern_name):
#     series_decomposition = decompose_series(dataframe[col_to_extract], id_patterns[pattern_name])
#     return series_decomposition


# def env_package_nomralizastion(dataframe, col_to_normalize, pattern_name):
#     dataframe[['lhs', 'rhs']] = dataframe[col_to_normalize].str.split('.', expand=True)
#     flag = dataframe['rhs'].apply(lambda x: x is None)
#     # antiflag = ~flag
#     temp = dataframe['lhs'][flag]
#     dataframe.loc[flag, 'rhs'] = temp
#     dataframe.loc[flag, 'lhs'] = ''
#     series_decomposition = decompose_series(dataframe['rhs'], id_patterns[pattern_name])
#     dataframe = pd.concat([dataframe, series_decomposition], axis=1)
#     # dataframe.append(series_decomposition, ignore_index=True)
#     # extracted = extract_with_pattern(dataframe, 'rhs', pattern_name)
#     return dataframe


def env_package_nomralizastion(dataframe, col_to_normalize, pattern_name, id_replacement_rule):
    dataframe[['lhs', 'rhs']] = dataframe[col_to_normalize].str.split('.', expand=True)
    flag = dataframe['rhs'].apply(lambda x: x is None)
    temp = dataframe['lhs'][flag]
    dataframe.loc[flag, 'rhs'] = temp
    dataframe.loc[flag, 'lhs'] = ''
    series_decomposition = decompose_series(dataframe['rhs'], id_replacement_rule)
    dataframe = pd.concat([dataframe, series_decomposition], axis=1)
    return dataframe


def add_overrides(dataframe, incol, outcol, override_dict):
    dataframe[outcol] = dataframe[incol]
    for key, value in override_dict.items():
        flag = dataframe[incol] == key
        dataframe.loc[flag, outcol] = value
    return dataframe


def flag_canonical(dataframe, incol, outcol, canonicals):
    dataframe[outcol] = False
    flag = dataframe[incol].isin(canonicals)
    dataframe.loc[flag, outcol] = True
    return dataframe


# ----


# # DATA
# biosample_sqlite_file = "/Users/MAM/Documents/gitrepos/biosample-analysis/target/harmonized_table.db"
# # TODO process these as a list?
# ncbitaxon_sqlite_file = "/Users/MAM/Documents/gitrepos/semantic-sql/db/ncbitaxon.db"
# envo_sqlite_file = "/Users/MAM/Documents/gitrepos/semantic-sql/db/envo.db"
# ncbitaxon_cnx = sqlite3.connect(ncbitaxon_sqlite_file)
# envo_cnx = sqlite3.connect(envo_sqlite_file)
# target_onto_prefix = 'ENVO'
# chars_to_whiteout = '._-'
# my_query_fields = ''
# my_row_req = 3
#
# env_package_overrides = {
#     'built environment': 'built',
#     'misc environment': 'miscellaneous',
#     'missing': 'no environmental package',
#     'unknown': 'no environmental package',
#     'default': 'no environmental package',
#     'unspecified': 'no environmental package',
#     'not available': 'no environmental package',
#     'not collected': 'no environmental package'
# }
#
# # Sample of the data we're working with
# biosample_cnx = sqlite3.connect(biosample_sqlite_file)
# q = """
# select
#     id,
#     env_package,
#     package,
#     package_name,
#     host_taxid,
#     taxonomy_id,
#     env_broad_scale,
#     env_local_scale,
#     env_medium
#     from biosample b
# limit 10
# """
# biosample_first_ten = pd.read_sql(q, biosample_cnx)
# print(biosample_first_ten)

# # Get the canonical checklist and package terms from NCBI
# # This document doesn't do a very good job of differentiating checklists (MIMAG, MIMARKS, etc.)
# # from packages (soil, water, etc.)
# package_dictionary = get_package_dictionary()
# print(package_dictionary)
# package_dictionary.to_sql('package_dictionary', biosample_cnx, if_exists='replace', index=False)

# # Do the Biosample checklist/package fields match any of the cannonical values?
# # How many Biosample rows are there?
# q = """
# select count(*) as biosample_row_count
# from biosample b
# """
# [biosample_row_count, query_duration] = timed_query(q, biosample_cnx, print_timing=False)
# print(biosample_row_count)
# print(query_duration)

# # How many of those rows can be inner-joined with the canonical checklists/packages?
# # Specifically, joining biosample.package_name = package_dictionary.DisplayName
# # TODO add indexing to docs and or makefile
# # create index biosample_package_name_idx on biosample(package_name);
# # create index package_dictionary_DisplayName_idx on package_dictionary(DisplayName);
# # create index biosample_package_idx on biosample(package);
# # create index biosample_p_pn_idx on biosample(package, package_name);
# q = """
# select
#     count(*) as cannonical_package_name_count
# from
#     biosample b
# inner join package_dictionary pd on
#     b.package_name = pd.DisplayName
# """
# [cannonical_package_name_count, query_duration] = timed_query(q, biosample_cnx, print_timing=True)
# print(cannonical_package_name_count)
# print(query_duration)

# # What do the combinations of package and package_name look like in the Biosample dataset?
# q = """
# select
#     package,
#     package_name,
#     count(*) as count
# from
#     biosample b
# group by
#     package ,
#     package_name
# order by
#     package ,
#     package_name
# """
# [package_name_combos, query_duration] = timed_query(q, biosample_cnx, print_timing=True)
# print(package_name_combos)
# print(query_duration)

# # What about the Biosample env_package values?
# # Are they also a small, highly regular set?
# q = """
# select
#     env_package,
#     count(*) as count
# from
#     biosample b
# group by
#     env_package
# order by
#     count(*) desc
# """
# [env_package_count, query_duration] = timed_query(q, biosample_cnx)
# print(env_package_count)
# print(query_duration)

# # env_package is going to need some cleanup
# # First, get a set of all canonical env_package values
# package_dictionary = make_tidy_col(package_dictionary, 'EnvPackage', 'eptidy')
# package_dictionary = make_tidy_col(package_dictionary, 'EnvPackageDisplay', 'epdtidy')
# # update in sqlite
# package_dictionary.to_sql('package_dictionary', biosample_cnx, if_exists='replace', index=False)
# valid_combo = []
# valid_combo = add_unique_to_list(valid_combo, package_dictionary['eptidy'])
# valid_combo = add_unique_to_list(valid_combo, package_dictionary['epdtidy'])
# print(valid_combo)

# # determine ID patterns
# q = """
# select
#     distinct stanza
#     from statements s
# where
#     predicate = 'rdf:type'
#     and "object" = 'owl:Class'
#     and stanza = subject"""
# # include non-envo IDs that come from envo?
# [ids_from_envo, query_duration] = timed_query(q, envo_cnx)
# ids_from_envo = add_prefix_col(ids_from_envo, 'stanza', 'prefix')
# id_patterns = get_multi_term_patterns(ids_from_envo, 'stanza', 'prefix')
# env_package_normalized = env_package_nomralizastion(env_package_count, 'env_package', target_onto_prefix)
# env_package_normalized = add_overrides(env_package_normalized, 'remaining_tidied', 'ep_override', env_package_overrides)
# env_package_normalized = flag_canonical(env_package_normalized, 'ep_override', 'is_canonical', valid_combo)
# env_package_normalized.to_sql('env_package_normalized', biosample_cnx, if_exists='replace', index=False)
# print(query_duration)

# ---- DONE-ISH?

# # What do the successful normalizations look like?
# q = """
# select
#     env_package,
#     count,
#     lhs,
#     extract,
#     ep_override
# from
#     env_package_normalized
# where
#     is_canonical = 1
# """
# [successful_normalizastions, query_duration] = timed_query(q, biosample_cnx)
# print(successful_normalizastions)
# print(query_duration)
#
# # Are there any normalization failures?
# q = """
# select
#     env_package,
#     count,
#     lhs,
#     extract,
#     ep_override
# from
#     env_package_normalized
# where
#     is_canonical = 0
# """
# [normalizastion_failures, query_duration] = timed_query(q, biosample_cnx)
# print(normalizastion_failures)
# print(query_duration)
#
# # utilizing ncbtitaxon
# # specifically, flag the biosamples whose taxon_id indicates they are an unclassified entity
# # ignoring the others will throw out samples OF multicellular organisms, like fruit flies
# # Add previous notes about what kinds of samples are missed by this bifurcation
# # like bacteria.unclassified_bacteria
#
# q = """
# select
#     distinct s.subject
# from
#     entailed_edge ee
# join statements s on
#     ee.subject = s.subject
# where
#     ee.predicate = 'rdfs:subClassOf'
#     and ee.object = 'NCBITaxon:2787823'
#     and s.predicate = 'rdfs:label'
# """
# [unclassified_taxa, query_duration] = timed_query(q, ncbitaxon_cnx)
# # unclassified_taxa['bare_id'] = unclassified_taxa['subject'].str.replace('NCBITaxon:', '', regex=True)
# # unclassified_taxa = unclassified_taxa[['bare_id']]
# # unclassified_taxa.columns = ['bare_id']
# unclassified_taxa['unclassified'] = True
# print(unclassified_taxa)
# print(query_duration)
#
# # SLOW
# # CHECK INDICES
# # Get the taxonomy_id values from the Biosamples
# q = """
# select
#     taxonomy_id biosample_taxid,
#     count(*) as count
# from
#     biosample b
# group by
#     taxonomy_id
# order by
#     count(*) desc
# """
# [biosample_tax_id_counts, query_duration] = timed_query(q, biosample_cnx)
# biosample_tax_id_counts['curie'] = 'NCBITaxon:' + biosample_tax_id_counts['biosample_taxid'].astype(str)
#
# # Merge the two taxon id datasets
# # I.e. flag the the Biosample records whose taxonomy_id field belongs to a subclass of 'unclassified entries'.
# biosample_tax_id_counts = biosample_tax_id_counts.merge(unclassified_taxa, left_on='curie',
#                                                         right_on='subject', how='left')
# biosample_tax_id_counts.unclassified.fillna(False, inplace=True)
# print(query_duration)
#
# # should really add labels to all of them
# q = """
# select
#     subject ,
#     value
# from statements
# where
#     predicate = 'rdfs:label' and subject = stanza
# """
# [all_tax_labels, query_duration] = timed_query(q, ncbitaxon_cnx)
#
# biosample_tax_id_counts = biosample_tax_id_counts.merge(all_tax_labels, left_on='curie',
#                                                         right_on='subject', how='left')
#
# biosample_tax_id_counts = biosample_tax_id_counts[['curie', 'biosample_taxid', 'count', 'unclassified', 'value']]
# biosample_tax_id_counts.columns = ['curie', 'biosample_taxid', 'count', 'unclassified', 'label']
# print(biosample_tax_id_counts)
# print(query_duration)
# biosample_tax_id_counts.to_sql('biobiosample_tax_id_counts', biosample_cnx, if_exists='replace', index=False)
#
# # Almost all of the taxa that are common in th biosample collection are either unclassified/metagenomes
# # or easily recognized cellular organisms
# # exceptions include:
# # 32630 = synthetic construct (other entries; other sequences; artificial sequences)
# #     'other entries' would add 16k rows on top of the 1k 'unclassified entities'
# #     metagenomes account for 331 of the 'unclassified entities'
# #     there are also a small number of uncultured/unclassified microorganisms in the biosample dataset
# # 77133 = uncultured bacterium (cellular organisms; Bacteria; environmental samples)
# #     'cellular organisms' would add 2M rows on top of the 1k 'unclassified entities'
# #     'cellular organisms; Bacteria; environmental samples' adds 26k
#
# # Get a table of scoped mixs annotations to be mapped to ontology classes.
# biosample_col_to_map = 'env_broad_scale'
# scoping_col = 'env_package_normalized.ep_override'
# scoping_value = 'water'
# # In this case, the scoping includes an inner join requirement for 'unclassified entities'
#
# q = 'select ' + biosample_col_to_map + """, count(*) as count
# from
#     biosample b
# join env_package_normalized on
#     b.env_package = env_package_normalized.env_package
# inner join biobiosample_tax_id_counts stic on
#     b.taxonomy_id = stic.biosample_taxid
# where """ + scoping_col + " = '" + scoping_value + \
#     "' group by " + biosample_col_to_map + """
# order by
#     count(*) desc"""
# [mapping_candidates, query_duration] = timed_query(q, biosample_cnx)
# print(mapping_candidates)
#
# # The Biosample format allows for pipe-delimited environmental package lists
# # Separate those out into their components
# multi_frames = []
# for row in mapping_candidates.itertuples(index=True, name='Pandas'):
#     split_check = row.env_broad_scale
#     if split_check is None:
#         split_check = ''
#     splitted = pd.Series(split_check.split("|"))
#     splitted_count = len(splitted)
#     repeated = [split_check] * splitted_count
#     repeated = pd.Series(repeated)
#     as_frame = pd.DataFrame(dict(repeated=repeated, splitted=splitted)).reset_index()
#     multi_frames.append(as_frame)
# concat_frame = pd.concat(multi_frames)
# concat_frame = concat_frame[['repeated', 'splitted']]
# mapping_candidates = mapping_candidates.merge(concat_frame, left_on='env_broad_scale', right_on='repeated', how='left')
# print(mapping_candidates)
#
# # Now try to extract ontology terms that are already present
# candidate_series_decomposition = decompose_series(mapping_candidates['splitted'], id_patterns[target_onto_prefix])
# mapping_candidates = pd.concat([mapping_candidates, candidate_series_decomposition], axis=1)
#
# # And join the extracted IDs with their labels
# ontodb = '/Users/MAM/Documents/gitrepos/semantic-sql/db/' + target_onto_prefix.lower() + '.db'
# ontocon = sqlite3.connect(ontodb)
# q = """
# select
#     subject ,
#     value
# from
#     statements s
# where
#     predicate = 'rdfs:label'
# """
# [onto_labels, query_duration] = timed_query(q, ontocon)
# mapping_candidates = mapping_candidates.merge(onto_labels, left_on='extract', right_on='subject', how='left')
#
# # Use cosine string distance to see if the labels match
# # I.e. the labels claimed by the Biosample data set and the labels asserted in the ontology
# # if they're close enough, consider the assigned ID legit
# # how close is close enough?
# my_cosine_obj = Cosine(1)
# mapping_candidates.value = mapping_candidates.value.fillna('')
# mapping_candidates['cosine'] = mapping_candidates.apply(
#     lambda my_row: my_cosine_obj.distance(my_row['remaining_tidied'].lower(), my_row['value'].lower()), axis=1)
#
# # Get ready to join in the other direction
# # I.e. trying to find ontology term IDs based on perfect label matches. Be careful not to reuse column names.
# mapping_candidates.columns = ['env_broad_scale', 'count', 'repeated', 'splitted', 'string', 'extract',
#                               'remaining_string', 'remaining_tidied', 'term_id', 'lab_from_id', 'lfi_cosine']
# mapping_candidates = mapping_candidates.merge(onto_labels, left_on='remaining_tidied', right_on='value', how='left')
# mapping_candidates.columns = ['env_broad_scale', 'count', 'repeated', 'splitted', 'string', 'extract',
#                               'remaining_string', 'remaining_tidied', 'term_id', 'lab_from_id',
#                               'lfi_cosine', 'term_id_from_lab', 'value']
#
# # Record a consensus
# # If either merging on codes or labels was successful.
# # cosines for first pass check on assigned IDs still haven't been filtered?
# mapping_candidates['consensus_id'] = mapping_candidates['term_id']
# mapping_candidates['consensus_id'][mapping_candidates['consensus_id'].isnull()] = \
#     mapping_candidates['term_id_from_lab'][mapping_candidates['consensus_id'].isnull()]
# mapping_candidates['consensus_lab'] = mapping_candidates['lab_from_id']
# mapping_candidates['consensus_lab'][mapping_candidates['consensus_lab'] == ''] = \
#     mapping_candidates['value'][mapping_candidates['consensus_lab'] == '']
# mapping_candidates.to_sql('mapping_scratch', biosample_cnx, if_exists='replace', index=False)
#
# # For which Biosample annotations were not mappings by merging found?
# # It looks like remaining_tidied is retaining too much punctuation
# # and loosing useful digits (relative to remaining_string)?
# # Should try harder to parse not-quite-right embedded IDs like ...
# needs_search = mapping_candidates.remaining_tidied[mapping_candidates.consensus_id.isna()]
# needs_search_counts = needs_search.value_counts()
# print(needs_search_counts)
#
# # ----
#
# # Use a search engine
# # For the mixs annotations that didn't already have cannonical IDs or labels
# ebs_raw_list = list(needs_search_counts.index)
# # get whiteout frame and relateds
# ebs_wo_frame = get_whiteout_frame(ebs_raw_list, replaced_chars=chars_to_whiteout)
# ebs_wo_list = get_wo_list(ebs_wo_frame)
# ebs_search_res = search_get_annotations_wrapper(ebs_wo_list, bad_chars=chars_to_whiteout, cat_name=biosample_col_to_map,
#                                                 ontoprefix=target_onto_prefix.lower(), query_fields='', rr=5)
# my_best_acceptable = get_best_acceptable(ebs_search_res)
# no_acceptable_mappings = get_no_acceptable_mappings(ebs_search_res, my_best_acceptable)
# # Some broad scales look like place names
# # Some get a good hit if 'biome' is added
# # how to manually review and then add back in?
# # add tp database
# #    no acceptable
# #    best acceptable
# #    ebs_search_results (no acceptable + all acceptable)?
# #    mapping_candidates -> mapping_scratch
