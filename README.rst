scoped-mapping
==============

When mapping input strings from column/field ``X`` in some datasource, consider the values in column/field ``Y`` to determine what ontology to map to.

Note that GitHub uses a hyphen and PyPI uses an underscore

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
  curl https://raw.githubusercontent.com/turbomam/scoped-mapping/main/requirements.txt -o requirements.txt
  curl https://raw.githubusercontent.com/turbomam/scoped-mapping/main/data/webmap_enums.yaml -o webmap_enums.yaml
  pip3 install -r requirements.txt
  pip3 install -i https://test.pypi.org/simple/ scoped-mapping


Sample code
-----------

Taking a list of strings as input::

  import scoped_mapping
  searchres_annotations = scoped_mapping.search_get_annotations_wrapper(
      ['Homo-sapiens', 'mus    musculus', 'rattus norvegicus'],
      bad_chars='._-', cat_name='test',  ontoprefix='ncbitaxon', query_fields='label')

Taking the enums from a LinkML model as input::

  import scoped_mapping
  import yaml
  import sys
  my_model_file = 'webmap_enums.yaml'
  my_selected_enum = 'Taxon_enum'
  my_model = scoped_mapping.read_yaml_model(my_model_file)
  yaml.safe_dump(my_model, sys.stdout, default_flow_style=False)
  yaml_mapped = scoped_mapping.map_from_yaml(my_model, my_selected_enum, print_enums=True, cat_name='unknown', ontoprefix='ncbitaxon')
  
Finding the best matches, reporting the faioures, and fedding the best matches back into the LinkML model::

  my_best_acceptable = scoped_mapping.get_best_acceptable(yaml_mapped)
  no_acceptable_mappings = scoped_mapping.get_no_acceptable_mappings(yaml_mapped, my_best_acceptable)
  scoped_mapping.rewrite_yaml(my_model, my_selected_enum, my_best_acceptable)
  yaml.safe_dump(my_model, sys.stdout, default_flow_style=False)
  


Scoping mappings based on subsets of NCBItaxon
----------------------------------------------

If a dataset has taxon values, one can use them to subset or scope how other values in the dataset should be mapped. For example, the NCBI Biosample metadata collection has MIxS triads (broad, narrow and medium) that could me mapped to ENVO terms in many cases. But ENVO might not be appropriate for cultured samples or samples that were taken from a multicellular organism. One way to check for those cases is looking for transitive subclasses in NCBItaxon. There are numerous ways to do that, but they are all generally computationally expensive.

Here, we use rdftab and relation-graph (via semantic-sql) to infer those transitive subClassOf relationships and load them into an SQLite database. Building this database requires lots of RAM and roughly 10 GB of disk space, but after that the querying is fast and convenient. Other ontologies can be loaded into the same database after that, but it might not be easy to infer thier transitive subclasses.

Our process of loading OWL ontologies into SQLite currently requires the `riot` library from Apache Jena. We don't have a suggested method for installing that into Windows or Linux yet. On MacOS, Jena can easily be installed with homebrew https://brew.sh/

``brew install jena``

::

  git clone git@github.com:cmungall/semantic-sql.git
  cd semantic-sql
  mkdir bin
  make bin/rdftab

Running subsequent steps may require editing semantic-sql's ``Makefile``. Specifically, invocation of relation-graph should be prepended with ``bin/``. For example::

  # ---
  # Inferences
  # ---
  # We use relation-graph
  inferences/%-inf.ttl: owl/%.owl
          bin/relation-graph 
        
etc.


``make bin/relation-graph``

Use as much RAM as possible. Otherwise inferring edges in NCBItaxon will be slow::

  export JAVA_OPTS=-Xmx24G
  make db/ncbitaxon.db
  
Now adding the ENVO ontology to that SQLite database is fast::

  make owl/envo.owl
  bin/rdftab db/ncbitaxon.db < owl/envo.owl 

Deploying
---------

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

