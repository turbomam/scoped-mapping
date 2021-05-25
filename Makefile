htdb_fileId = 1r_FOSBBNa5qXAd3Upu2V--VdivlVehL3
htdb_fileName = target/harmonized_table.db.gz
/tmp/htdb_fileId:
	# get code for confirming download of large file without virus scanning
	# TODO switch to using temporary files instead of named files that somebody else might delete or overwrite
	[ ! -e /tmp/gdcookie ] || rm /tmp/gdcookie
	[ ! -e /tmp/htdb_fileId ] || rm /tmp/htdb_fileId
	curl -sc /tmp/gdcookie "https://drive.google.com/uc?export=download&id=$(htdb_fileId)" > /dev/null; \
	awk '/_warning_/ {print $$NF}' /tmp/gdcookie > /tmp/htdb_fileId

harmonized_table.db: /tmp/htdb_fileId
# gets a gzipped SQLite database
# that was derived from biosample metadata records
# that have a "harmonized name" field
# in ftp://ftp.ncbi.nlm.nih.gov/biosample/biosample_set.xml.gz
# via Google Drive
	$(eval gc4m=$(shell cat /tmp/htdb_fileId))
	curl -Lb /tmp/gdcookie "https://drive.google.com/uc?export=download&confirm=$(gc4m)&id=$(htdb_fileId)" -o $(htdb_fileName)
	# is force really a good idea here?
	# or just make another rule?
	# does Bill already have one?
	ls -lh target/harmonized_table.db*
	gunzip -f $(htdb_fileName)
	ls -lh target/harmonized_table.db*
	# test for existence approach
	[ ! -e /tmp/gdcookie ] || rm /tmp/gdcookie
	[ ! -e /tmp/htdb_fileId ] || rm /tmp/htdb_fileId
	sqlite3 target/harmonized_table.db < ht_indicies.sql

semantic-sql:
	# requires Apache Jena's riot library
	# following MacOS homebrew installation approach
	git clone git@github.com:cmungall/semantic-sql.git
	cd semantic-sql ; mkdir bin ; \
		sed -i .bak 's^relation-graph --ontology-file^bin/relation-graph --ontology-file^' Makefile ; \
		make bin/rdftab ; \
		make bin/relation-graph; \
		brew install jena

semantic-sql/db/ncbitaxon.db: semantic-sql
	date
	export JAVA_OPTS=-Xmx24G ; \
	cd semantic-sql ; \
	make db/ncbitaxon.db
	date

semantic-sql/db/envo.db: semantic-sql
	date
	export JAVA_OPTS=-Xmx24G ; \
	cd semantic-sql ; \
	make db/envo.db
	date

.PHONY all: harmonized_table.db semantic-sql/db/ncbitaxon.db semantic-sql/db/envo.db

clean:
	rm -rf harmonized_table.db
	rm -rf semantic-sql

