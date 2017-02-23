"""
EMBL format (:mod:`skbio.io.format.embl`)
===============================================

.. currentmodule:: skbio.io.format.embl

EMBL format (EMBL Flat File Format) stores sequence and its
annotation together. Each line begins with a two-character line code, 
which indicates the type of information contained in the line. The currently 
used line types, along with their respective line codes, are listed below:

    +--------------------------------+----------------------------------+
    | ID - identification            | (begins each entry; 1 per entry) |
    +--------------------------------+----------------------------------+
    | AC - accession number          | (>=1 per entry)                  |
    +--------------------------------+----------------------------------+
    | PR - project identifier        | (0 or 1 per entry)               |
    +--------------------------------+----------------------------------+
    | DT - date                      | (2 per entry)                    |
    +--------------------------------+----------------------------------+
    | DE - description               | (>=1 per entry)                  |
    +--------------------------------+----------------------------------+
    | KW - keyword                   | (>=1 per entry)                  |
    +--------------------------------+----------------------------------+
    | OS - organism species          | (>=1 per entry)                  |
    +--------------------------------+----------------------------------+
    | OC - organism classification   | (>=1 per entry)                  |
    +--------------------------------+----------------------------------+
    | OG - organelle                 | (0 or 1 per entry)               |
    +--------------------------------+----------------------------------+
    | RN - reference number          | (>=1 per entry)                  |
    +--------------------------------+----------------------------------+
    | RC - reference comment         | (>=0 per entry)                  |
    +--------------------------------+----------------------------------+
    | RP - reference positions       | (>=1 per entry)                  |
    +--------------------------------+----------------------------------+
    | RX - reference cross-reference | (>=0 per entry)                  |
    +--------------------------------+----------------------------------+
    | RG - reference group           | (>=0 per entry)                  |
    +--------------------------------+----------------------------------+
    | RA - reference author(s)       | (>=0 per entry)                  |
    +--------------------------------+----------------------------------+
    | RT - reference title           | (>=1 per entry)                  |
    +--------------------------------+----------------------------------+
    | RL - reference location        | (>=1 per entry)                  |
    +--------------------------------+----------------------------------+
    | DR - database cross-reference  | (>=0 per entry)                  |
    +--------------------------------+----------------------------------+
    | CC - comments or notes         | (>=0 per entry)                  |
    +--------------------------------+----------------------------------+
    | AH - assembly header           | (0 or 1 per entry)               |
    +--------------------------------+----------------------------------+
    | AS - assembly information      | (0 or >=1 per entry)             |
    +--------------------------------+----------------------------------+
    | FH - feature table header      | (2 per entry)                    |
    +--------------------------------+----------------------------------+
    | FT - feature table data        | (>=2 per entry)                  |
    +--------------------------------+----------------------------------+
    | XX - spacer line               | (many per entry)                 |
    +--------------------------------+----------------------------------+
    | SQ - sequence header           | (1 per entry)                    |
    +--------------------------------+----------------------------------+
    | CO - contig/construct line     | (0 or >=1 per entry)             |
    +--------------------------------+----------------------------------+
    | bb - (blanks) sequence data    | (>=1 per entry)                  |
    +--------------------------------+----------------------------------+
    | // - termination line          | (ends each entry; 1 per entry)   |
    +--------------------------------+----------------------------------+

Note that some entries will not contain all of the line types, and some line
types occur many times in a single entry. As indicated, each entry begins with
an identification line (ID) and ends with a terminator line (//). The various 
line types appear in entries in the order in which they are listed above 
(except for XX lines which may appear anywhere between the ID and SQ lines). 
An example of a EMBL file can be seen here [1]_.

In general, fixed format items have been kept to a 
minimum, and a more syntax-oriented structure adopted for the lines. 
The two exceptions to this are the sequence data lines and the feature table
lines, for which a fixed format was felt to offer significant advantages
to the user. Users who write programs to process the database entries should
not make any assumptions about the column placement of items on lines other
than these two: all other line types are free-format. 

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.sequence.Sequence`                                 |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.DNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.RNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.Protein`                                  |
+------+------+---------------------------------------------------------------+
|Yes   | Yes  | generator of :mod:`skbio.sequence.Sequence` objects           |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
**State: Experimental as of 0.6.**

Sections before Feature Table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
All the sections before Feature Table will be read into the attribute
of ``metadata``. The header and its content of a section is stored as
a pair of key and value in ``metadata``. For the Reference
section, its value is stored as a list, as there are often multiple
reference sections in one EMBL record.

Feature Table section
^^^^^^^^^^^^^^^^^^^^
The International Nucleotide Sequence Database Collaboration (INSDC
[2]_) is a joint effort among the DDBJ, EMBL-EBI, and NCBI. These
organisations all use the same "Feature Table" layout in their plain
text flat file formats, which are documented in detail [3]_. The
feature keys and their qualifiers are also described in this webpage
[4]_.

The Feature Table section will be stored in ``interval_metadata`` of
``Sequence`` or its sub-class. Each sub-section is stored as an
``Interval`` object in ``interval_metadata``. Each ``Interval`` object
has ``metadata`` keeping the information of this feature in the
sub-section.

To normalize the vocabulary between multiple formats (currently only
the INSDC Feature Table and GFF3) to store metadata of interval
features, we rename some terms in some formats to the same common name
when parsing them into memory, as described in this table:

+-----------+-----------+-----------+---------+------------------------------+
|INSDC      |GFF3       |Key stored |Value    |Description                   |
|feature    |columns or |           |type     |                              |
|table      |attributes |           |stored   |                              |
+===========+===========+===========+=========+==============================+
|inference  |source     |source     |str      |the algorithm or experiment   |
|           |(column 2) |           |         |used to generate this feature |
+-----------+-----------+-----------+---------+------------------------------+
|feature key|type       |type       |str      |the type of the feature       |
|           |(column 3) |           |         |                              |
+-----------+-----------+-----------+---------+------------------------------+
|N/A        |score      |score      |float    |the score of the feature      |
|           |(column 6) |           |         |                              |
+-----------+-----------+-----------+---------+------------------------------+
|N/A        |strand     |strand     |str      |the strand of the feature. +  |
|           |(column 7) |           |         |for positive strand, - for    |
|           |           |           |         |minus strand, and . for       |
|           |           |           |         |features that are not         |
|           |           |           |         |stranded. In addition, ?  can |
|           |           |           |         |be used for features whose    |
|           |           |           |         |strandedness is relevant, but |
|           |           |           |         |unknown.                      |
+-----------+-----------+-----------+---------+------------------------------+
|codon_start|phase      |phase      |int      |the offset at which the first |
|           |(column 8) |           |         |complete codon of a coding    |
|           |           |           |         |feature can be found, relative|
|           |           |           |         |to the first base of that     |
|           |           |           |         |feature. It is 0, 1, or 2 in  |
|           |           |           |         |GFF3 or 1, 2, or 3 in EMBL.|
|           |           |           |         |The stored value is 0, 1, or  |
|           |           |           |         |2, following in GFF3 format.  |
+-----------+-----------+-----------+---------+------------------------------+
|db_xref    |Dbxref     |db_xref    |list of  |A database cross reference    |
|           |           |           |str      |                              |
+-----------+-----------+-----------+---------+------------------------------+
|N/A        |ID         |ID         |str      |feature ID                    |
+-----------+-----------+-----------+---------+------------------------------+
|note       |Note       |note       |str      |any comment or additional     |
|           |           |           |         |information                   |
+-----------+-----------+-----------+---------+------------------------------+
|translation|N/A        |translation|str      |the protein sequence for CDS  |
|           |           |           |         |features                      |
+-----------+-----------+-----------+---------+------------------------------+

``Location`` string
+++++++++++++++++++
There are 5 types of location descriptors defined in Feature
Table. This explains how they will be parsed into the bounds of
``Interval`` object (note it converts the 1-based coordinate to
0-based):

    1. a single base number. e.g. 67. This is parsed to ``(66, 67)``.

    2. a site between two neighboring bases. e.g. 67^68. This
       is parsed to ``(66, 67)``.

    3. a single base from inside a range. e.g. 67.89. This is parsed to
       ``(66, 89)``.

    4. a pair of base numbers defining a sequence span. e.g. 67..89. This
       is parsed to ``(66, 89)``.

    5. a remote sequence identifier followed by a location descriptor
       defined above. e.g. J00123.1:67..89. This will be discarded
       because it is not on the current sequence. When it is combined
       with local descriptor like J00123.1:67..89,200..209, the
       local part will be kept to be ``(199, 209)``.


.. note:: The Location string is fully stored in ``Interval.metadata``
   with key ``__location``.  The key starting with ``__`` is "private"
   and should be modified with care.


Sequence section
^^^^^^^^^^^^^^^^^^
The sequence in the Sequence section is always in lowercase for
the EMBL files downloaded from EBI. For the RNA molecules, ``t``
(thymine), instead of ``u`` (uracil) is used in the sequence. All
EMBL writers follow these conventions while writing EMBL files.


Format Parameters
-----------------

Reader-specific Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^
The ``constructor`` parameter can be used with the ``Sequence`` generator
to specify the in-memory type of each EMBL record that is parsed.
``constructor`` should be ``Sequence`` or a sub-class of ``Sequence``.
It is also detected by the unit label on the LOCUS line. For example, if it
is ``bp``, it will be read into ``DNA``; if it is ``aa``, it will be read
into ``Protein``. Otherwise, it will be read into ``Sequence``. This default
behavior is overridden by setting ``constructor``.

``lowercase`` is another parameter available for all EMBL readers.
By default, it is set to ``True`` to read in the sequence
as lowercase letters. This parameter is passed to ``Sequence`` or
its sub-class constructor.

``seq_num`` is a parameter used with the ``Sequence``, ``DNA``, ``RNA``, and
``Protein`` EMBL readers. It specifies which EMBL record to read from
a EMBL file with multiple records in it.

Examples
--------

Reading and Writing EMBL Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Suppose we have the following EMBL file example modified from [5]_:

>>> eb_str = '''
... ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.
... XX
... AC   X56734; S46826;
... XX
... DT   12-SEP-1991 (Rel. 29, Created)
... DT   25-NOV-2005 (Rel. 85, Last updated, Version 11)
... XX
... DE   Trifolium repens mRNA for non-cyanogenic beta-glucosidase
... XX
... KW   beta-glucosidase.
... XX
... OS   Trifolium repens (white clover)
... OC   Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
... OC   Spermatophyta; Magnoliophyta; eudicotyledons; Gunneridae; Pentapetalae;
... OC   rosids; fabids; Fabales; Fabaceae; Papilionoideae; Trifolieae; Trifolium.
... XX
... RN   [5]
... RP   1-1859
... RX   DOI; 10.1007/BF00039495.
... RX   PUBMED; 1907511.
... RA   Oxtoby E., Dunn M.A., Pancoro A., Hughes M.A.;
... RT   "Nucleotide and derived amino acid sequence of the cyanogenic
... RT   beta-glucosidase (linamarase) from white clover (Trifolium repens L.)";
... RL   Plant Mol. Biol. 17(2):209-219(1991).
... XX
... RN   [6]
... RP   1-1859
... RA   Hughes M.A.;
... RT   ;
... RL   Submitted (19-NOV-1990) to the INSDC.
... RL   Hughes M.A., University of Newcastle Upon Tyne, Medical School, Newcastle
... RL   Upon Tyne, NE2 4HH, UK
... XX
... DR   MD5; 1e51ca3a5450c43524b9185c236cc5cc.
... XX
... FH   Key             Location/Qualifiers
... FH
... FT   source          1..1859
... FT                   /organism="Trifolium repens"
... FT                   /mol_type="mRNA"
... FT                   /clone_lib="lambda gt10"
... FT                   /clone="TRE361"
... FT                   /tissue_type="leaves"
... FT                   /db_xref="taxon:3899"
... FT   mRNA            1..1859
... FT                   /experiment="experimental evidence, no additional details
... FT                   recorded"
... FT   CDS             14..1495
... FT                   /product="beta-glucosidase"
... FT                   /EC_number="3.2.1.21"
... FT                   /note="non-cyanogenic"
... FT                   /db_xref="GOA:P26204"
... FT                   /db_xref="InterPro:IPR001360"
... FT                   /db_xref="InterPro:IPR013781"
... FT                   /db_xref="InterPro:IPR017853"
... FT                   /db_xref="InterPro:IPR033132"
... FT                   /db_xref="UniProtKB/Swiss-Prot:P26204"
... FT                   /protein_id="CAA40058.1"
... FT                   /translation="MDFIVAIFALFVISSFTITSTNAVEASTLLDIGNLSRSSFPRGFI
... FT                   FGAGSSAYQFEGAVNEGGRGPSIWDTFTHKYPEKIRDGSNADITVDQYHRYKEDVGIMK
... FT                   DQNMDSYRFSISWPRILPKGKLSGGINHEGIKYYNNLINELLANGIQPFVTLFHWDLPQ
... FT                   VLEDEYGGFLNSGVINDFRDYTDLCFKEFGDRVRYWSTLNEPWVFSNSGYALGTNAPGR
... FT                   CSASNVAKPGDSGTGPYIVTHNQILAHAEAVHVYKTKYQAYQKGKIGITLVSNWLMPLD
... FT                   DNSIPDIKAAERSLDFQFGLFMEQLTTGDYSKSMRRIVKNRLPKFSKFESSLVNGSFDF
... FT                   IGINYYSSSYISNAPSHGNAKPSYSTNPMTNISFEKHGIPLGPRAASIWIYVYPYMFIQ
... FT                   EDFEIFCYILKINITILQFSITENGMNEFNDATLPVEEALLNTYRIDYYYRHLYYIRSA
... FT                   IRAGSNVKGFYAWSFLDCNEWFAGFTVRFGLNFVD"
... XX
... SQ   Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
...      aaacaaacca aatatggatt ttattgtagc catatttgct ctgtttgtta ttagctcatt        60
...      cacaattact tccacaaatg cagttgaagc ttctactctt cttgacatag gtaacctgag       120
...      tcggagcagt tttcctcgtg gcttcatctt tggtgctgga tcttcagcat accaatttga       180
...      aggtgcagta aacgaaggcg gtagaggacc aagtatttgg gataccttca cccataaata       240
...      tccagaaaaa ataagggatg gaagcaatgc agacatcacg gttgaccaat atcaccgcta       300
...      caaggaagat gttgggatta tgaaggatca aaatatggat tcgtatagat tctcaatctc       360
...      ttggccaaga atactcccaa agggaaagtt gagcggaggc ataaatcacg aaggaatcaa       420
...      atattacaac aaccttatca acgaactatt ggctaacggt atacaaccat ttgtaactct       480
...      ttttcattgg gatcttcccc aagtcttaga agatgagtat ggtggtttct taaactccgg       540
...      tgtaataaat gattttcgag actatacgga tctttgcttc aaggaatttg gagatagagt       600
...      gaggtattgg agtactctaa atgagccatg ggtgtttagc aattctggat atgcactagg       660
...      aacaaatgca ccaggtcgat gttcggcctc caacgtggcc aagcctggtg attctggaac       720
...      aggaccttat atagttacac acaatcaaat tcttgctcat gcagaagctg tacatgtgta       780
...      taagactaaa taccaggcat atcaaaaggg aaagataggc ataacgttgg tatctaactg       840
...      gttaatgcca cttgatgata atagcatacc agatataaag gctgccgaga gatcacttga       900
...      cttccaattt ggattgttta tggaacaatt aacaacagga gattattcta agagcatgcg       960
...      gcgtatagtt aaaaaccgat tacctaagtt ctcaaaattc gaatcaagcc tagtgaatgg      1020
...      ttcatttgat tttattggta taaactatta ctcttctagt tatattagca atgccccttc      1080
...      acatggcaat gccaaaccca gttactcaac aaatcctatg accaatattt catttgaaaa      1140
...      acatgggata cccttaggtc caagggctgc ttcaatttgg atatatgttt atccatatat      1200
...      gtttatccaa gaggacttcg agatcttttg ttacatatta aaaataaata taacaatcct      1260
...      gcaattttca atcactgaaa atggtatgaa tgaattcaac gatgcaacac ttccagtaga      1320
...      agaagctctt ttgaatactt acagaattga ttactattac cgtcacttat actacattcg      1380
...      ttctgcaatc agggctggct caaatgtgaa gggtttttac gcatggtcat ttttggactg      1440
...      taatgaatgg tttgcaggct ttactgttcg ttttggatta aactttgtag attagaaaga      1500
...      tggattaaaa aggtacccta agctttctgc ccaatggtac aagaactttc tcaaaagaaa      1560
...      ctagctagta ttattaaaag aactttgtag tagattacag tacatcgttt gaagttgagt      1620
...      tggtgcacct aattaaataa aagaggttac tcttaacata tttttaggcc attcgttgtg      1680
...      aagttgttag gctgttattt ctattatact atgttgtagt aataagtgca ttgttgtacc      1740
...      agaagctatg atcataacta taggttgatc cttcatgtat cagtttgatg ttgagaatac      1800
...      tttgaattaa aagtcttttt ttattttttt aaaaaaaaaa aaaaaaaaaa aaaaaaaaa       1859
... //
... '''


Now we can read it as ``DNA`` object:

>>> import io
>>> from skbio import DNA, RNA, Sequence
>>> eb = io.StringIO(eb_str)
>>> dna_seq = DNA.read(eb)
>>> dna_seq
DNA
----------------------------------------------------------------------
Metadata:
    'AC': <class 'list'>
    'DE': <class 'list'>
    'DR': <class 'list'>
    'DT': <class 'list'>
    'ID': <class 'dict'>
    'KW': <class 'list'>
    'OC': <class 'list'>
    'OS': <class 'list'>
    'RN': <class 'list'>
    'SQ': <class 'dict'>
Interval metadata:
    3 interval features
Stats:
    length: 1859
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 35.99%
----------------------------------------------------------------------
0    AAACAAACCA AATATGGATT TTATTGTAGC CATATTTGCT CTGTTTGTTA TTAGCTCATT
60   CACAATTACT TCCACAAATG CAGTTGAAGC TTCTACTCTT CTTGACATAG GTAACCTGAG
...
1740 AGAAGCTATG ATCATAACTA TAGGTTGATC CTTCATGTAT CAGTTTGATG TTGAGAATAC
1800 TTTGAATTAA AAGTCTTTTT TTATTTTTTT AAAAAAAAAA AAAAAAAAAA AAAAAAAAA


Since this is a riboswitch molecule, we may want to read it as
``RNA``.  As the EMBL file usually have ``t`` instead of ``u`` in
the sequence, we can read it as ``RNA`` by converting ``t`` to ``u``:

>>> eb = io.StringIO(eb_str)
>>> rna_seq = RNA.read(eb)
>>> rna_seq
RNA
----------------------------------------------------------------------
Metadata:
    'AC': <class 'list'>
    'DE': <class 'list'>
    'DR': <class 'list'>
    'DT': <class 'list'>
    'ID': <class 'dict'>
    'KW': <class 'list'>
    'OC': <class 'list'>
    'OS': <class 'list'>
    'RN': <class 'list'>
    'SQ': <class 'dict'>
Interval metadata:
    3 interval features
Stats:
    length: 1859
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 35.99%
----------------------------------------------------------------------
0    AAACAAACCA AAUAUGGAUU UUAUUGUAGC CAUAUUUGCU CUGUUUGUUA UUAGCUCAUU
60   CACAAUUACU UCCACAAAUG CAGUUGAAGC UUCUACUCUU CUUGACAUAG GUAACCUGAG
...
1740 AGAAGCUAUG AUCAUAACUA UAGGUUGAUC CUUCAUGUAU CAGUUUGAUG UUGAGAAUAC
1800 UUUGAAUUAA AAGUCUUUUU UUAUUUUUUU AAAAAAAAAA AAAAAAAAAA AAAAAAAAA

>>> rna_seq == dna_seq.transcribe()
True

References
----------
.. [1] http://www.ebi.ac.uk/ena/data/view/BN000065&display=text
.. [2] http://www.insdc.org/
.. [3] http://www.insdc.org/files/feature_table.html
.. [4] http://www.ebi.ac.uk/ena/WebFeat/
.. [5] http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?X56734

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re
from functools import partial
from itertools import chain

from skbio.io import create_format, EMBLFormatError
from skbio.io.format._base import (
    _get_nth_sequence, _line_generator, _too_many_blanks)
from skbio.util._misc import chunk_str
from skbio.sequence import Sequence, DNA, RNA, Protein
from skbio.io.format._sequence_feature_vocabulary import (
    _yield_section, _parse_section_default, _serialize_section_default,
    _parse_feature_table, _serialize_feature_table)


embl = create_format('embl')

# This list is ordered. From EMBL specification
_HEADERS = [('ID', 'ID'),  # identification            (begins each entry; 1 per entry)
            ('AC', 'AC'),  # accession number          (>=1 per entry)
            ('PR', 'PR'),  # project identifier        (0 or 1 per entry)
            ('DT', 'DT'),  # date                      (2 per entry)
            ('DE', 'DR'),  # description               (>=1 per entry)
            ('KW', 'KW'),  # keyword                   (>=1 per entry)
            ('OS', 'OS'),  # organism species          (>=1 per entry)
            ('OC', 'OC'),  # organism classification   (>=1 per entry)
            ('OG', 'OG'),  # organelle                 (0 or 1 per entry)
            ('RN', 'RN'),  # reference number          (>=1 per entry)
            ('RC', 'RN'),  # reference comment         (>=0 per entry)
            ('RP', 'RN'),  # reference positions       (>=1 per entry)
            ('RX', 'RN'),  # reference cross-reference (>=0 per entry)
            ('RG', 'RN'),  # reference group           (>=0 per entry)
            ('RA', 'RN'),  # reference author(s)       (>=0 per entry)
            ('RT', 'RN'),  # reference title           (>=1 per entry)
            ('RL', 'RN'),  # reference location        (>=1 per entry)
            ('DR', 'DR'),  # database cross-reference  (>=0 per entry)
            ('CC', 'CC'),  # comments or notes         (>=0 per entry)
            ('AS', 'AS'),  # assembly information      (0 or >=1 per entry)
            ('FT', 'FT'),  # feature table data        (>=2 per entry)
            ('SQ', 'SQ'),  # sequence header           (1 per entry)
            ('CO', 'CO')]  # contig/construct line     (0 or >=1 per entry)

@embl.sniffer()
def _embl_sniffer(fh):
    # check the 1st real line is a valid ID line
    #print('embl.sniffer')
    if _too_many_blanks(fh, 5):
        return False, {}
    try:
        line = next(_line_generator(fh, skip_blanks=True, strip=False))
    except StopIteration:
        return False, {}

    # TODO: differentiate between Swiss-Prot and EMBL
    #       EMBL: ID   <1>; SV <2>; <3>; <4>; <5>; <6>; <7> BP.
    # Swiss-Prot: ID   EntryName Status; SequenceLength.
    if line.startswith('ID'):
        return True, {}
    else:
        return False, {}


@embl.reader(None)
def _embl_to_generator(fh, constructor=None, **kwargs):
    for record in _parse_records(fh, _parse_single_embl):
        yield _construct(record, constructor, **kwargs)


@embl.reader(Sequence)
def _embl_to_sequence(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh, _parse_single_embl), seq_num)
    return _construct(record, Sequence, **kwargs)


@embl.reader(DNA)
def _embl_to_dna(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh, _parse_single_embl), seq_num)
    return _construct(record, DNA, **kwargs)


@embl.reader(RNA)
def _embl_to_rna(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh, _parse_single_embl), seq_num)
    return _construct(record, RNA, **kwargs)


@embl.reader(Protein)
def _embl_to_protein(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh, _parse_single_embl), seq_num)
    return _construct(record, Protein, **kwargs)


# @embl.writer(None)
# def _generator_to_embl(obj, fh):
#     for obj_i in obj:
#         _serialize_single_embl(obj_i, fh)


# @embl.writer(Sequence)
# def _sequence_to_embl(obj, fh):
#     _serialize_single_embl(obj, fh)


# @embl.writer(DNA)
# def _dna_to_embl(obj, fh):
#     _serialize_single_embl(obj, fh)


# @embl.writer(RNA)
# def _rna_to_embl(obj, fh):
#     _serialize_single_embl(obj, fh)


# @embl.writer(Protein)
# def _protein_to_embl(obj, fh):
#     _serialize_single_embl(obj, fh)


def _construct(record, constructor=None, **kwargs):
    '''Construct the object of Sequence, DNA, RNA, or Protein.
    '''
    seq, md, imd = record
    if 'lowercase' not in kwargs:
        kwargs['lowercase'] = True
    if constructor is None:
        unit = md['LOCUS']['unit']
        if unit == 'bp':
            constructor = DNA
        elif unit == 'aa':
            constructor = Protein

    if constructor == RNA:
        return DNA(
            seq, metadata=md, interval_metadata=imd, **kwargs).transcribe()
    else:
        return constructor(
            seq, metadata=md, interval_metadata=imd, **kwargs)


def _parse_records(fh, parser):
    data_chunks = []
    for line in _line_generator(fh, skip_blanks=True, strip=False):
        if line.startswith('//'):
            yield parser(data_chunks)
            data_chunks = []
        else:
            data_chunks.append(line)


# dictionary to determine which headers belong in the same section
_HEADER_SECTION = dict(_HEADERS)
def _header_section(header):
    return _HEADER_SECTION.get(header, header)

def _parse_single_embl(chunks):
    metadata = {}
    interval_metadata = None
    sequence = ''

    # each section starts with a different HEADER.
    section_splitter = _yield_section(
        lambda x, i, l: (not x[0].isspace()) and
        (i > 0 and (_header_section(x[:2]) != _header_section(l[i - 1][:2]))),
        strip=False)

    for section in section_splitter(chunks):
        header = section[0].split(None, 1)[0]
        if header not in _HEADER_SECTION: # ignore XX, AH, FH
            continue
        parser = _PARSER_TABLE.get(header, _parse_section_default)

        if header == 'FT':
            # This requires 'ID' line parsed before 'FT', which should
            # be true and is implicitly checked by the sniffer.
            parsed = parser(
                [l[2:] for l in section],
                length=metadata['ID']['size'],
                feature_indent=' ' * 19)
        else:
            parsed = parser(section)

        # reference can appear multiple times
        if header == 'RN':
            if header in metadata:
                metadata[header].append(parsed)
            else:
                metadata[header] = [parsed]
        elif header == 'SQ':
            metadata[header] = parsed[0]
            sequence = parsed[1]
        elif header == 'FT':
            interval_metadata = parsed
        else:
            metadata[header] = parsed
    return sequence, metadata, interval_metadata


# def _serialize_single_embl(obj, fh):
#     '''Write a EMBL record.

#     Always write it in EBI canonical way:
#     1. sequence in lowercase
#     2. 'u' as 't' even in RNA molecules.

#     Parameters
#     ----------
#     obj : Sequence or its child class

#     '''
#     # write out the headers
#     md = obj.metadata
#     for header in _HEADERS:
#         serializer = _SERIALIZER_TABLE.get(
#             header, _serialize_section_default)
#         if header in md:
#             out = serializer(header, md[header])
#             # test if 'out' is a iterator.
#             # cf. Effective Python Item 17
#             if iter(out) is iter(out):
#                 for s in out:
#                     fh.write(s)
#             else:
#                 fh.write(out)
#         if header == 'FEATURES':
#             if obj.has_interval_metadata():
#                 # magic number 21: the amount of indentation before
#                 # feature table starts as defined by INSDC
#                 indent = 21
#                 fh.write('{header:<{indent}}Location/Qualifiers\n'.format(
#                     header=header, indent=indent))
#                 for s in serializer(obj.interval_metadata._intervals, indent):
#                     fh.write(s)
#     # write out the sequence
#     # always write RNA seq as DNA
#     if isinstance(obj, RNA):
#         obj = obj.reverse_transcribe()

#     # always write in lowercase
#     seq_str = str(obj).lower()

#     for s in _serialize_origin(seq_str):
#         fh.write(s)
#     fh.write('//\n')

def _parse_id(lines):
    '''Parse ID line.
    
    Example:
    ID   CD789012; SV 4; linear; genomic DNA; HTG; MAM; 500 BP.
    ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.
    
    The ID (IDentification) line is always the first line of an entry. The
    format of the ID line is:
    ID   <1>; SV <2>; <3>; <4>; <5>; <6>; <7> BP.
    The tokens represent:
    1. Primary accession number
    2. Sequence version number
    3. Topology: "circular", "linear"
    4. Molecule type: "genomic DNA", "genomic RNA", "mRNA", "tRNA", "rRNA", 
       "other RNA", "other DNA", "transcribed RNA", "viral cRNA", 
       "unassigned DNA", "unassigned RNA"
    5. Data class: "CON", "PAT", "EST", "GSS", "HTC", "HTG", "MGA", 
       "WGS", "TSA", "STS", "STD"
    6. Taxonomic division: "PHG", "ENV", "FUN", "HUM", "INV", "MAM", 
       "VRT", "MUS", "PLN", "PRO", "ROD", "SYN", "TGN", "UNC", "VRL"
    7. Sequence length [1]_

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt

    Test
    ----
    line = ('ID   CD789012; SV 4; linear; genomic DNA; HTG; MAM; 500 BP.\n'
    'ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.')
    '''
    pattern = (r'ID   '
               r'(?P<accession>[^\s;]+); '
               r'SV (?P<seq_version>[^\s;]+); '
               r'(?P<topology>linear|circular); '
               r'(?P<mol_type>genomic DNA|genomic RNA|mRNA|tRNA|rRNA|'
               r'other RNA|other DNA|transcribed RNA|viral cRNA|'
               r'unassigned DNA|unassigned RNA); '
               r'(?P<data_class>CON|PAT|EST|GSS|HTC|HTG|MGA|WGS|'
               r'TSA|STS|STD); '
               r'(?P<tax_division>PHG|ENV|FUN|HUM|INV|MAM|VRT|MUS|'
               r'PLN|PRO|ROD|SYN|TGN|UNC|VRL); '
               r'(?P<size>\d+) BP.')
    line = lines[0]
    matches = re.match(pattern, line)
    if not matches:
        raise EMBLFormatError(
            "Could not parse the ID line:\n%s" % line)
    res = matches.groupdict()
    res['size'] = int(res['size'])
    return res

def _parse_ac(lines):
    '''Parse AC line. (>=1 per entry)
    The AC (ACcession number) line lists the accession numbers associated with 
    the entry.
    Each accession number, or range of accession numbers, is terminated by a
    semicolon. Where necessary, more than one AC line is used. Consecutive
    secondary accession numbers in ENA flatfiles are shown in the form of 
    inclusive accession number ranges [1]_.

    Example:
    AC   X56734; S46826;
    AC   Y00001; X00001-X00005; X00008; Z00001-Z00005;

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    return list(
        chain(*[[ac.rstrip(';') for ac in line[5:].split()]
                for line in lines]))


def _parse_pr(lines):
    '''Parse PR line. (0 or 1 per entry)
    The PR (PRoject) line shows the International Nucleotide Sequence Database
    Collaboration (INSDC) Project Identifier that has been assigned to the entry.
    Full details of INSDC Project are available at
    http://www.ebi.ac.uk/ena/about/page.php?page=project_guidelines [1]_.

    Example:
    PR   Project:17285;

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    if len(lines) > 1:
        raise EMBLFormatError(
            "Expected 0 or 1 PR lines, found {:d}.".format(len(line)))
    return lines[0][5:].strip().rstrip(';')

def _parse_dt(lines):
    '''Parse DT line. (2 per entry)
    The DT (DaTe) line shows when an entry first appeared in the database and
    when it was last updated.  Each entry contains two DT lines, formatted
    as follows:
    DT   DD-MON-YYYY (Rel. #, Created)
    DT   DD-MON-YYYY (Rel. #, Last updated, Version #)
    [1]_.

    Example:
    DT   12-SEP-1991 (Rel. 29, Created)
    DT   13-SEP-1993 (Rel. 37, Last updated, Version 8)

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    if len(lines) != 2:
        raise EMBLFormatError(
            "Expected 2 DT lines, found {:d}.".format(len(lines)))
    return [line[5:] for line in lines]

def _parse_de(lines):
    '''Parse DE line. (>=1 per entry)
    The description is given in ordinary English and is free-format. Often, more
    than one DE line is required; when this is the case, the text is divided only
    between words.
    The first DE line generally contains a brief description, which can stand
    alone for cataloguing purposes. [1]_.

    Example:
    DE   Trifolium repens mRNA for non-cyanogenic beta-glucosidase

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    return [line[5:] for line in lines]

def _parse_kw(lines):
    '''Parse KW line. (>=1 per entry)
    The KW (KeyWord) lines provide information which can be used to generate
    cross-reference indexes of the sequence entries based on functional,
    structural, or other categories deemed important.
    The format for a KW line is:
        KW   keyword[; keyword ...].
    More than one keyword may be listed on each KW line; the keywords are
    separated by semicolons, and the last keyword is followed by a full
    stop. Keywords may consist of more than one word, and they may contain
    embedded blanks and stops. A keyword is never split between lines.
    The keywords are ordered alphabetically; the ordering implies no hierarchy
    of importance or function.  If an entry has no keywords assigned to it,
    it will contain a single KW line like this:
        KW   .
    [1]_.

    Example:
    KW   beta-glucosidase.

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    return list(
        chain(*[[ac.strip().rstrip('.') for ac in line[5:].split(';')]
                for line in lines]))

def _parse_os(lines):
    '''Parse OS line. (>=1 per entry)
    The OS (Organism Species) line specifies the preferred scientific name of
    the organism which was the source of the stored sequence. In most 
    cases this is done by giving the Latin genus and species designations, 
    followed (in parentheses) by the preferred common name in English where
    known. The format is:
        OS   Genus species (name)
    [1]_.

    Example:
    OS   Trifolium repens (white clover)
    OS   Mus musculus x Rattus norvegicus

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    return [line[5:] for line in lines]

def _parse_oc(lines):
    '''Parse OC line. (>=1 per entry)
    The OC (Organism Classification) lines contain the taxonomic classification
    Of the source organism. 
    The classification is listed top-down as nodes in a taxonomic tree in which 
    the most general grouping is given first.  The classification may be 
    distributed over several OC lines, but nodes are not split or hyphenated 
    between lines. The individual items are separated by semicolons and the
    list is terminated by a full stop. The format for the OC line is:
        OC   Node[; Node...].
    [1]_.

    Example:
    OC   Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
    OC   euphyllophytes; Spermatophyta; Magnoliophyta; eudicotyledons; Rosidae;
    OC   Fabales; Fabaceae; Papilionoideae; Trifolium.

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    return list(
        chain(*[[ac.strip().rstrip('.') for ac in line[5:].split(';')]
                for line in lines]))

def _parse_og(lines):
    '''Parse OG line. (0 or 1 per entry)
    The OG (OrGanelle) linetype indicates the sub-cellular location of non-nuclear
    sequences.  It is only present in entries containing non-nuclear sequences
    and appears after the last OC line in such entries.
    [1]_.

    Example:
    OS   Euglena gracilis (green algae)
    OC   Eukaryota; Planta; Phycophyta; Euglenophyceae.
    OG   Plastid:Chloroplast

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    if len(lines) > 1:
        raise EMBLFormatError(
            "Expected 0 or 1 OG lines, found {:d}.".format(len(line)))
    return lines[0][5:].strip()

def _parse_rn(lines):
    '''The Reference (RN, RC, RP, RX, RG, RA, RT, RL) Lines
    These lines comprise the literature citations within the database.
    The citations provide access to the papers from which the data has been 
    abstracted. The reference lines for a given citation occur in a block, and
    are always in the order RN, RC, RP, RX, RG, RA, RT, RL. Within each such 
    reference block the RN line occurs once, the RC, RP and RX lines occur zero
    or more times, and the RA, RT, RL lines each occur one or more times. 
    If several references are given, there will be a reference block for each. [1]_.

    Example:
    RN   [2]
    RP   1-1657990
    RG   Prochlorococcus genome consortium
    RA   Larimer F.;
    RT   ;
    RL   Submitted (03-JUL-2003) to the INSDC.
    RL   Larimer F., DOE Joint Genome Institute, Production Genomics Facility, 
    RL   2800 Mitchell Drive, Walnut Creek, CA 94598, USA, and the Genome 
    RL   Analysis Group, Oak Ridge National Laboratory, 1060 Commerce Park Drive, 
    RL   Oak Ridge, TN 37831, USA;

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    res = {
        'RN': lines[0][5:],
        'RC': [],
        'RP': [],
        'RX': [],
        'RG': [],
        'RA': [],
        'RT': [],
        'RL': []
    }
    for line in lines[1:]:
        res[line[:2]].append(line[5:])

    return res

def _parse_dr(lines):
    '''Parse DR line. (>=0 per entry)
    The DR (Database Cross-reference) line cross-references other databases which
    contain information related to the entry in which the DR line appears.
    Format:
    DR   database_identifier; primary_identifier; secondary_identifier.
    [1]_.

    Example:
    DR   MGI; 98599; Tcrb-V4.

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    return [[t.rstrip(';.') for t in line[5:].split()] for line in lines]

def _parse_cc(lines):
    '''Parse CC line. (>=0 per entry)
    CC lines are free text comments about the entry, and may be used to convey 
    any sort of information thought to be useful that is unsuitable for
    inclusion in other line types.
    [1]_.

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    return [line[5:] for line in lines]

def _parse_as(lines):
    '''Parse AS line. (0 or >=1 per entry)
    The AS (ASsembly Information) lines provide information on the composition of 
    a TPA or TSA sequence. These lines include information on local sequence spans
    (those spans seen in the sequence of the entry showing the AS lines) plus
    identifiers and base spans of contributing primary sequences (for ENA
    primary entries only).
    [1]_.

    Example:
    AH   LOCAL_SPAN     PRIMARY_IDENTIFIER     PRIMARY_SPAN     COMP
    AS   1-426          AC004528.1             18665-19090         
    AS   427-526        AC001234.2             1-100            c
    AS   527-1000       TI55475028             not_available

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    return [line[5:].split() for line in lines]

def _parse_sq(lines):
    '''Parse SQ line. (1 per entry)
    The SQ (SeQuence header) line marks the beginning of the sequence data and 
    Gives a summary of its content.

    The sequence data line has a line code consisting of two blanks. The sequence
    is written 60 bases per line, in groups of 10 bases separated by a blank
    character, beginning at position 6 of the line. The direction listed is 
    always 5' to 3', and wherever possible the non-coding strand 
    (homologous to the message) has been stored. Columns 73-80 of each 
    sequence line contain base numbers for easier reading and quick 
    location of regions of interest. The numbers are right justified and indicate
    the number of the last base on each line.
    [1]_.

    Example:
    SQ   Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
        aaacaaacca aatatggatt ttattgtagc catatttgct ctgtttgtta ttagctcatt        60
        cacaattact tccacaaatg cagttgaagc ttctactctt cttgacatag gtaacctgag       120

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    pattern = (r'SQ   '
               r'Sequence (?P<total>\d+) BP; '
               r'(?P<a>\d+) A; '
               r'(?P<c>\d+) C; '
               r'(?P<g>\d+) G; '
               r'(?P<t>\d+) T; '
               r'(?P<other>\d+) other;')
    line = lines[0]
    matches = re.match(pattern, line)
    if not matches:
        raise EMBLFormatError(
            "Could not parse the SQ line:\n%s" % line)
    summary = dict([(k, int(v)) for k, v in matches.groupdict().items()])
    # remove the base number at the end of each data line
    sequence = ''.join(chain(*[line.split()[:-1] for line in lines[1:]]))
    return (summary, sequence)

def _parse_co(lines):
    '''Parse CO line. (0 or >=1 per entry)
    Con(structed) sequences in the CON data classes represent complete
    chromosomes, genomes and other long sequences constructed from segment entries.
    CON data class entries do not contain sequence data per se, but rather the
    assembly information on all accession.versions and sequence locations relevant
    to building the constructed sequence. The assembly information is represented in
    the CO lines.
    [1]_.

    Example:
    CO   join(Z99104.1:1..213080,Z99105.1:18431..221160,Z99106.1:13061..209100, 
    CO   Z99107.1:11151..213190,Z99108.1:11071..208430,Z99109.1:11751..210440, 
    CO   Z99110.1:15551..216750,Z99111.1:16351..208230,Z99112.1:4601..208780, 
    CO   Z99113.1:26001..233780,Z99114.1:14811..207730,Z99115.1:12361..213680, 
    CO   Z99116.1:13961..218470,Z99117.1:14281..213420,Z99118.1:17741..218410, 
    CO   Z99119.1:15771..215640,Z99120.1:16411..217420,Z99121.1:14871..209510, 
    CO   Z99122.1:11971..212610,Z99123.1:11301..212150,Z99124.1:11271..215534) 

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
    '''
    return [line[5:] for line in lines]



_PARSER_TABLE = {
    'ID': _parse_id,
    'AC': _parse_ac,
    'PR': _parse_pr,
    'DT': _parse_dt,
    'DE': _parse_de,
    'KW': _parse_kw,
    'OS': _parse_os,
    'OC': _parse_oc,
    'OG': _parse_og,
    'RN': _parse_rn, # includes RN, RC, RP, RX, RG, RA, RT, RL
    'DR': _parse_dr,
    'CC': _parse_cc,
    'AS': _parse_as,
    'FT': _parse_feature_table,
    'SQ': _parse_sq,
    'CO': _parse_co
}

