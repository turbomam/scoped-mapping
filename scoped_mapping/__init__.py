
import sys
import requests as requests
from pkg_resources import get_distribution, DistributionNotFound
import pandas as pd
import re
import urllib
from strsimpy.cosine import Cosine
import yaml

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    pass  # package is not installed

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


def adder(x, y):
    return x + y


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
    # pp.pprint(userdata)
    resp = ols_session.get('https://www.ebi.ac.uk/ols/api/search', params=userdata)
    rj = resp.json()
    # pp.pprint(list(rj.keys()))
    # pp.pprint(rj['responseHeader'])
    rjf = pd.DataFrame(rj['response']['docs'])
    rjf_rows = len(rjf.index)
    if rjf_rows < 1:
        col_count = len(ols_cols)
        blank_row = pd.Series([''] * col_count)
        rjf = rjf.append(blank_row, ignore_index=True)
        rjf.columns = ols_cols
    # clogger.info(rjf)
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
        search_annotations_merge.apply(lambda row:
                                       cosine_obj.distance(row.query.strip().lower(),
                                                           row['name'].strip().lower()), axis=1)
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
                                   ontoprefix='', query_fields=''):
    my_category_search_results_frame = get_csr_frame(raw_list, bad_chars, category_name=cat_name,
                                                     ontoprefix=ontoprefix, query_fields=query_fields)
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
    for row in best_acceptable.itertuples(index=True, name='Pandas'):
        model['enums'][enum]['permissible_values'][row.raw]['meaning'] = row.obo_id
        model['enums'][enum]['permissible_values'][row.raw]['description'] = row.label


# my_model_file = '/Users/MAM/webmap_enums.yaml'
# my_selected_enum = 'Taxon_enum'
# my_model = read_yaml_model(my_model_file)
# yaml_mapped = map_from_yaml(my_model, my_selected_enum, print_enums=True, cat_name='unknown', ontoprefix='ncbitaxon')
# my_best_acceptable = get_best_acceptable(yaml_mapped)
# no_acceptable_mappings = get_no_acceptable_mappings(yaml_mapped, my_best_acceptable)
# rewrite_yaml(my_model, my_selected_enum, my_best_acceptable)
# yaml.safe_dump(my_model, sys.stdout, default_flow_style=False)
