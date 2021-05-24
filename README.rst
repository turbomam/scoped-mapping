scoped-mapping
==============

When mapping input strings from column/field ``X`` in some datasource, consider the values in column/field ``Y`` to determine what ontology to map to.

Note that GitHub uses a hyphen and PyPI uses an underscore

Currently compatible with MacOS. Requires the riot library from Apache Jena. ``make all`` uses homebrew for installing Jena, but does not install homebrew.



Installation
------------
::

  python3.9 -m venv sm_venv
  source sm_venv/bin/activate
  pip3 install -r requirements.txt
  pip3 install -i https://test.pypi.org/simple/ scoped-mapping
  make all
  


Sample code
-----------

See Jupyter Notebooks  


Scoping mappings based on subsets of NCBItaxon
----------------------------------------------

If a dataset has taxon values, one can use them to subset or scope how other values in the dataset should be mapped. For example, the NCBI Biosample metadata collection has MIxS triads (broad, narrow and medium) that could me mapped to ENVO terms in many cases. But ENVO might not be appropriate for cultured samples or samples that were taken from a multicellular organism. One way to check for those cases is looking for transitive subclasses in NCBItaxon. There are numerous ways to do that, but they are all generally computationally expensive.

Here, we use rdftab and relation-graph (via semantic-sql) to infer those transitive subClassOf relationships and load them into an SQLite database. Building this database requires lots of RAM and roughly 10 GB of disk space, but after that the querying is fast and convenient.

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

