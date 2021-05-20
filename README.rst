scoped-mapping
==============

When mapping input strings from column/field ``X`` in some datasource, consider the values in column/field ``Y`` to determine what ontology to map to.

Note that GitHub uses a hypen and PyPI uses an underscore

Dependencies
------------
- Python 3.9?
- pip3?
- venv?
- setuptools>=56.2.0
- pandas>=1.2.4
- requests>=2.25.1
- strsimpy>=0.2.0
- PyYAML>=5.4.1




Installation
------------
::

  python3.9 -m venv sm_venv
  source sm_venv/bin/activate
  pip3 install pandas PyYAML requests strsimpy
  pip3 install -i https://test.pypi.org/simple/ scoped-mapping


Sample code
-----------

taking a list of strings as input::

  import scoped_mapping
  searchres_annotations = scoped_mapping.search_get_annotations_wrapper(
      ['Homo-sapiens', 'mus    musculus', 'rattus norvegicus'],
      bad_chars='._-', cat_name='test',  ontoprefix='ncbitaxon', query_fields='label')

taking the enums from a LinkML model as input::

  import scoped_mapping
  my_model_file = 'data/webmap_enums.yaml'
  my_selected_enum = 'Taxon_enum'

  yaml_mapped = scoped_mapping.map_from_yaml(my_model_file, my_selected_enum, print_enums=True,
                                             cat_name='unknown', ontoprefix='ncbitaxon')
  my_best_acceptable = scoped_mapping.get_best_acceptable(yaml_mapped)
  no_acceptable_mappings = scoped_mapping.get_no_acceptable_mappings(yaml_mapped, my_best_acceptable)
  
----

Deploying

::

  git add ...
  git commit -m ...
  git push
  git tag ...
  python3.9 -m build --sdist --wheel .
  
rm dist/scoped_mapping... (old versions)

::

  twine upload --repository pypitest dist/*

uninstall and reinstall locally (see above)

