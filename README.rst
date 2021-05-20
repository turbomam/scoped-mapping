scoped-mapping
==============

When mapping input strings from column/field ``X`` in some datasource, consider the values in column/field ``Y`` to determine what ontology to map to.

Note that GitHub uses a hypen and PyPI uses an underscore

Dependencies
------------

Installation
------------

``pip3 install -i https://test.pypi.org/simple/ scoped-mapping``

Sample code
-----------

taking a list of strings as input::

  searchres_annotations = scoped_mapping.search_get_annotations_wrapper(
      ['Homo-sapiens', 'mus    musculus', 'rattus norvegicus'],
      bad_chars='._-', cat_name='test',  ontoprefix='ncbitaxon', query_fields='label')

taking the enums from a LinkML model as input::

  my_model_file = 'data/webmap_enums.yaml'
  my_selected_enum = 'Taxon_enum'

  yaml_mapped = scoped_mapping.map_from_yaml(my_model_file, my_selected_enum, print_enums=True,
                                             cat_name='unknown', ontoprefix='ncbitaxon')
  my_best_acceptable = scoped_mapping.get_best_acceptable(yaml_mapped)
  no_acceptable_mappings = scoped_mapping.get_no_acceptable_mappings(yaml_mapped, my_best_acceptable)
