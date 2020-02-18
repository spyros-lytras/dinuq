=====
DINUQ
=====


Usage
=====


Locally importing the package
-----------------------------

In the root directory do:
``pip install .``

Then in python do:
``import dinuq``


dinuq.SDU()
-----------

The SDU module will calculate the synonymous dinucleotide usage for all sequences in a given fasta file

*Arguments*

``sdu = dinuq.SDU(fasta_file, dinucl, position = ['bridge'], boots = 'none')``

``fasta file #required``

``dinucl = ['CpC', 'CpG', 'CpU', 'CpA', 'GpC', 'GpG', 'GpU', 'GpA', 'UpC', 'UpG', 'UpU', 'UpA', 'ApC', 'ApG', 'ApU', 'ApA'] #required``

``position = ['pos1', 'pos2', 'bridge'] #default is bridge``

``boots = integer #default is none``

*Output*

``sdu = {'accession': {'dinucleotideposition': [sdu_value, [bootstrap_value1, bootstrap_value2, bootstrap_valuen]]}}``


dinuq.RSDU()
------------

The RSDU module will calculate the relative synonymous dinucleotide usage for all sequences in a given fasta file

*Arguments*

``rsdu = dinuq.RSDU(fasta_file, dinucl, position = ['bridge'], boots = 'none')``

``fasta file #required``

``dinucl = ['CpC', 'CpG', 'CpU', 'CpA', 'GpC', 'GpG', 'GpU', 'GpA', 'UpC', 'UpG', 'UpU', 'UpA', 'ApC', 'ApG', 'ApU', 'ApA'] #required``

``position = ['pos1', 'pos2', 'bridge'] #default is bridge``

``boots = integer #default is none``

*Output*

``rsdu = {'accession': {'dinucleotideposition': [rsdu_value, [bootstrap_value1, bootstrap_value2, bootstrap_valuen]]}}``


dinuq.dict_to_tsv()
-------------------

This module creates a tsv of a table with the sdu or rsdu dictionaries, the user can choose what error to include (STDEV, SEM, MIN-MAX)

*Arguments*

``dinuq.dict_to_tsv(dictionary, output_file, error = 'none')``

``dictionary = sdu or rsdu #required``

``output_file #required``

``error = 'none', #default``

``	'extrema' #minimum and maximum of bootstrapped distribution``

``	'stdev' #mean plus/minus the distribution's standard deviation``

``	'sem' #mean plus/minus the distribution's standard error of the mean``
	
	
dinuq.RDA()
-----------

The RDA module will calculate the relative dinucleotide abundance for all sequences in a given fasta file, either for the entire sequence or specific positions

*Arguments*

``rda = dinuq.RDA(fasta_file, dinucl, position = ['all'])``

``fasta_file #required``

``dinucl = ['CpC', 'CpG', 'CpU', 'CpA', 'GpC', 'GpG', 'GpU', 'GpA', 'UpC', 'UpG', 'UpU', 'UpA', 'ApC', 'ApG', 'ApU', 'ApA'] #required``

``position = ['pos1', 'pos2', 'bridge', 'all'] #default is all``

*Output*

``rda = {'accession': {'dinucleotideposition': [rda_value]}}``	


dinuq.RDA_to_tsv()
-------------------

This module creates a tsv of a table with the rda dictionary

*Arguments*

``dinuq.RDA_to_tsv(dictionary, output_file)``

``dictionary = rda #required``

``output_file #required``


dinuq.RSCU()
-----------

The RSCU module will calculate the relative synonymous codon usage for all sequences in a given fasta file

*Arguments*

``rscu = dinuq.RSCU(fasta_file)``

``fasta_file #required``


*Output*

``rscu = {'accession': {'codon': rscu_value}}``


dinuq.RSCU_to_tsv()
-------------------

This module creates a tsv of a table with the rscu dictionary

*Arguments*

``dinuq.RSCU_to_tsv(dictionary, output_file)``

``dictionary = rscu #required``

``output_file #required``


