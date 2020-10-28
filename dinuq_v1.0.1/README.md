DinuQ
=====

The DinuQ (Dinucleotide Quantification) Python3 package provides a range
of metrics for quantifying dinucleotide representation and synonymous
codon usage in a DNA/RNA sequence. These include the recently developed
Synonymous Dinucleotide Usage (SDU) and Relative Synonymous Dinucleotide
Usage (RSDU).

Citation: *[Lytras, S.; Hughes, J. Synonymous Dinucleotide Usage: A Codon-Aware Metric for Quantifying Dinucleotide Representation in Viruses. Viruses 2020, 12, 462](https://doi.org/10.3390/v12040462)*

Usage
-----

### Package installation

Using pip, in a Unix terminal do: `pip install dinuq`

Then in python do: `import dinuq`

### Modules

#### dinuq.SDU()

The SDU module will calculate the Synonymous Dinucleotide Usage for all
sequences in a given fasta file.



*Arguments*

Required arguments:

-   a fasta file with:
    :   -   any number of coding sequences (no internal stop codons)
        -   a different, preferably short, fasta header for each
            sequence (e.g. an accession)

-   A list of dinucleotides of interest (still needs to be a list if
    it's only one, e.g. ['CpG'])

Optional arguments:

-   A list of dinucleotide frame positions. By default the module will
    only calculate the SDU for the bridge position, for each specified
    dinucleotide.
-   If you want to calculate error intervals for the SDU values, you can
    specify a number of iterations for the error measuring method
    (suggested value is 1000). Notice that this will significantly slow
    down the calculation.

`sdu = dinuq.SDU(fasta_file, dinucl, position = ['bridge'], samples = 'none')`

-   `fasta file #required`
-   `dinucl = ['CpC', 'CpG', 'CpU', 'CpA', 'GpC', 'GpG', 'GpU', 'GpA', 'UpC', 'UpG', 'UpU', 'UpA', 'ApC', 'ApG', 'ApU', 'ApA'] #required`
-   `position = ['pos1', 'pos2', 'bridge'] #default is bridge`
-   `samples = integer #default is none`



*Output*

The output of the module is a dictionary of accessions as keys and inner
dictionaries as values. The inner dictionaries have each dinucleotide
position as keys (e.g. CpGbridge) and a list of calculated SDU values as
the value. If the error margins are being calculated, an inner list of
SDU values calculated for each random sampling (specified in the samples
argument) is included.

`sdu = {'accession': {'dinucleotideposition': [sdu_value, [bootstrap_value1, bootstrap_value2, bootstrap_valuen]]}}`



#### dinuq.RSDU()

The RSDU module will calculate the Relative Synonymous Dinucleotide
Usage for all sequences in a given fasta file.



*Arguments*

The arguments are the same as the these for the SDU module.

`rsdu = dinuq.RSDU(fasta_file, dinucl, position = ['bridge'], samples = 'none')`

-   `fasta file #required`
-   `dinucl = ['CpC', 'CpG', 'CpU', 'CpA', 'GpC', 'GpG', 'GpU', 'GpA', 'UpC', 'UpG', 'UpU', 'UpA', 'ApC', 'ApG', 'ApU', 'ApA'] #required`
-   `position = ['pos1', 'pos2', 'bridge'] #default is bridge`
-   `samples = integer #default is none`



*Output*

The output format is the same as in the SDU module.

`rsdu = {'accession': {'dinucleotideposition': [rsdu_value, [bootstrap_value1, bootstrap_value2, bootstrap_valuen]]}}`



#### dinuq.dict\_to\_tsv()

This module creates a tsv file in your working directory with the sdu or
rsdu dictionary information in a table format. The user can choose how
to summarise the error distribution (STDEV, SEM, MIN-MAX) if that has
been calculated.



*Arguments*

Required arguments:

-   a sdu or rsdu dictionary produced by the SDU or RSDU module
    respectively
-   A name for the output tsv file

Optional arguments:

-   A summary of the error distribution (given that it has been calculated by the SDU/RSDU module). This can be:
    :   -   The minimum and maximum value of the distribution (extrema)
        -   The standard deviation margins around the error
            distribution's mean (stdev)
        -   The standard error of the mean margins around the mean (sem)

`dinuq.dict_to_tsv(dictionary, output_file, error = 'none')`

-   `dictionary = sdu or rsdu #required`
-   `output_file #required`
-   `error = 'none', #default`
    :   -   `'extrema' #minimum and maximum of bootstrapped distribution`
        -   `'stdev' #mean plus/minus the distribution's standard deviation`
        -   `'sem' #mean plus/minus the distribution's standard error of the mean`



#### dinuq.RDA()

The RDA module will calculate the Relative Dinucleotide Abundance for
all sequences in a given fasta file, either for the entire sequence or
specific dinucleotide frame positions.



*Arguments*

Required arguments:

-   a fasta file with:
    :   -   any number of coding sequences (no internal stop codons)
        -   a different, preferably short, fasta header for each
            sequence (e.g. an accession)

-   A list of dinucleotides of interest (still needs to be a list if
    it's only one, e.g. ['CpG'])

Optional arguments:

-   A list of dinucleotide frame positions. By default the module will
    calculate the RDA for the entire sequence (no frame position
    separation).

`rda = dinuq.RDA(fasta_file, dinucl, position = ['all'])`

-   `fasta_file #required`
-   `dinucl = ['CpC', 'CpG', 'CpU', 'CpA', 'GpC', 'GpG', 'GpU', 'GpA', 'UpC', 'UpG', 'UpU', 'UpA', 'ApC', 'ApG', 'ApU', 'ApA'] #required`
-   `position = ['pos1', 'pos2', 'bridge', 'all'] #default is all`



*Output*

The output of the module is a dictionary of accessions as keys and inner
dictionaries as values. The inner dictionaries have each dinucleotide
position as keys (e.g. CpGbridge) and a list of the calculated RDA value
as the value.

`rda = {'accession': {'dinucleotideposition': [rda_value]}}`



#### dinuq.RDA\_to\_tsv()

This module creates a tsv file in your working directory with the rda
dictionary information in a table format.



*Arguments*

Required arguments:

-   a rda dictionary produced by the RDA module
-   A name for the output tsv file

`dinuq.RDA_to_tsv(dictionary, output_file)`

`dictionary = rda #required`

`output_file #required`



#### dinuq.RSCU()

The RSCU module will calculate the Relative Synonymous Codon Usage for
all sequences in a given fasta file.



*Arguments*

Required arguments:

-   a fasta file with:
    :   -   any number of coding sequences (no internal stop codons)
        -   a different, preferably short, fasta header for each
            sequence (e.g. an accession)

`rscu = dinuq.RSCU(fasta_file)`

-   `fasta_file #required`



*Output*

The output of the module is a dictionary of accessions as keys and inner
dictionaries as values. The inner dictionaries have each codon as keys
and the calculated RSCU value as the value.

`rscu = {'accession': {'codon': rscu_value}}`



#### dinuq.RSCU\_to\_tsv()

This module creates a tsv file in your working directory with the rscu
dictionary information in a table format.



*Arguments*

Required arguments:

-   a rscu dictionary produced by the RSCU module
-   A name for the output tsv file

`dinuq.RSCU_to_tsv(dictionary, output_file)`

`dictionary = rscu #required`

`output_file #required`
