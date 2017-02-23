# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

from skbio import Protein, DNA, RNA, Sequence
from skbio.metadata import IntervalMetadata
from skbio.util import get_data_path
from skbio.io import EMBLFormatError
from skbio.io.format.embl import (
    _embl_sniffer)


eb_str = '''
ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.
XX
AC   X56734; S46826;
XX
DT   12-SEP-1991 (Rel. 29, Created)
DT   25-NOV-2005 (Rel. 85, Last updated, Version 11)
XX
DE   Trifolium repens mRNA for non-cyanogenic beta-glucosidase
XX
KW   beta-glucosidase.
XX
OS   Trifolium repens (white clover)
OC   Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
OC   Spermatophyta; Magnoliophyta; eudicotyledons; Gunneridae; Pentapetalae;
OC   rosids; fabids; Fabales; Fabaceae; Papilionoideae; Trifolieae; Trifolium.
XX
RN   [5]
RP   1-1859
RX   DOI; 10.1007/BF00039495.
RX   PUBMED; 1907511.
RA   Oxtoby E., Dunn M.A., Pancoro A., Hughes M.A.;
RT   "Nucleotide and derived amino acid sequence of the cyanogenic
RT   beta-glucosidase (linamarase) from white clover (Trifolium repens L.)";
RL   Plant Mol. Biol. 17(2):209-219(1991).
XX
RN   [6]
RP   1-1859
RA   Hughes M.A.;
RT   ;
RL   Submitted (19-NOV-1990) to the INSDC.
RL   Hughes M.A., University of Newcastle Upon Tyne, Medical School, Newcastle
RL   Upon Tyne, NE2 4HH, UK
XX
DR   MD5; 1e51ca3a5450c43524b9185c236cc5cc.
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..1859
FT                   /organism="Trifolium repens"
FT                   /mol_type="mRNA"
FT                   /clone_lib="lambda gt10"
FT                   /clone="TRE361"
FT                   /tissue_type="leaves"
FT                   /db_xref="taxon:3899"
FT   mRNA            1..1859
FT                   /experiment="experimental evidence, no additional details
FT                   recorded"
FT   CDS             14..1495
FT                   /product="beta-glucosidase"
FT                   /EC_number="3.2.1.21"
FT                   /note="non-cyanogenic"
FT                   /db_xref="GOA:P26204"
FT                   /db_xref="InterPro:IPR001360"
FT                   /db_xref="InterPro:IPR013781"
FT                   /db_xref="InterPro:IPR017853"
FT                   /db_xref="InterPro:IPR033132"
FT                   /db_xref="UniProtKB/Swiss-Prot:P26204"
FT                   /protein_id="CAA40058.1"
FT                   /translation="MDFIVAIFALFVISSFTITSTNAVEASTLLDIGNLSRSSFPRGFI
FT                   FGAGSSAYQFEGAVNEGGRGPSIWDTFTHKYPEKIRDGSNADITVDQYHRYKEDVGIMK
FT                   DQNMDSYRFSISWPRILPKGKLSGGINHEGIKYYNNLINELLANGIQPFVTLFHWDLPQ
FT                   VLEDEYGGFLNSGVINDFRDYTDLCFKEFGDRVRYWSTLNEPWVFSNSGYALGTNAPGR
FT                   CSASNVAKPGDSGTGPYIVTHNQILAHAEAVHVYKTKYQAYQKGKIGITLVSNWLMPLD
FT                   DNSIPDIKAAERSLDFQFGLFMEQLTTGDYSKSMRRIVKNRLPKFSKFESSLVNGSFDF
FT                   IGINYYSSSYISNAPSHGNAKPSYSTNPMTNISFEKHGIPLGPRAASIWIYVYPYMFIQ
FT                   EDFEIFCYILKINITILQFSITENGMNEFNDATLPVEEALLNTYRIDYYYRHLYYIRSA
FT                   IRAGSNVKGFYAWSFLDCNEWFAGFTVRFGLNFVD"
XX
SQ   Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
     aaacaaacca aatatggatt ttattgtagc catatttgct ctgtttgtta ttagctcatt        60
     cacaattact tccacaaatg cagttgaagc ttctactctt cttgacatag gtaacctgag       120
     tcggagcagt tttcctcgtg gcttcatctt tggtgctgga tcttcagcat accaatttga       180
     aggtgcagta aacgaaggcg gtagaggacc aagtatttgg gataccttca cccataaata       240
     tccagaaaaa ataagggatg gaagcaatgc agacatcacg gttgaccaat atcaccgcta       300
     caaggaagat gttgggatta tgaaggatca aaatatggat tcgtatagat tctcaatctc       360
     ttggccaaga atactcccaa agggaaagtt gagcggaggc ataaatcacg aaggaatcaa       420
     atattacaac aaccttatca acgaactatt ggctaacggt atacaaccat ttgtaactct       480
     ttttcattgg gatcttcccc aagtcttaga agatgagtat ggtggtttct taaactccgg       540
     tgtaataaat gattttcgag actatacgga tctttgcttc aaggaatttg gagatagagt       600
     gaggtattgg agtactctaa atgagccatg ggtgtttagc aattctggat atgcactagg       660
     aacaaatgca ccaggtcgat gttcggcctc caacgtggcc aagcctggtg attctggaac       720
     aggaccttat atagttacac acaatcaaat tcttgctcat gcagaagctg tacatgtgta       780
     taagactaaa taccaggcat atcaaaaggg aaagataggc ataacgttgg tatctaactg       840
     gttaatgcca cttgatgata atagcatacc agatataaag gctgccgaga gatcacttga       900
     cttccaattt ggattgttta tggaacaatt aacaacagga gattattcta agagcatgcg       960
     gcgtatagtt aaaaaccgat tacctaagtt ctcaaaattc gaatcaagcc tagtgaatgg      1020
     ttcatttgat tttattggta taaactatta ctcttctagt tatattagca atgccccttc      1080
     acatggcaat gccaaaccca gttactcaac aaatcctatg accaatattt catttgaaaa      1140
     acatgggata cccttaggtc caagggctgc ttcaatttgg atatatgttt atccatatat      1200
     gtttatccaa gaggacttcg agatcttttg ttacatatta aaaataaata taacaatcct      1260
     gcaattttca atcactgaaa atggtatgaa tgaattcaac gatgcaacac ttccagtaga      1320
     agaagctctt ttgaatactt acagaattga ttactattac cgtcacttat actacattcg      1380
     ttctgcaatc agggctggct caaatgtgaa gggtttttac gcatggtcat ttttggactg      1440
     taatgaatgg tttgcaggct ttactgttcg ttttggatta aactttgtag attagaaaga      1500
     tggattaaaa aggtacccta agctttctgc ccaatggtac aagaactttc tcaaaagaaa      1560
     ctagctagta ttattaaaag aactttgtag tagattacag tacatcgttt gaagttgagt      1620
     tggtgcacct aattaaataa aagaggttac tcttaacata tttttaggcc attcgttgtg      1680
     aagttgttag gctgttattt ctattatact atgttgtagt aataagtgca ttgttgtacc      1740
     agaagctatg atcataacta taggttgatc cttcatgtat cagtttgatg ttgagaatac      1800
     tttgaattaa aagtcttttt ttattttttt aaaaaaaaaa aaaaaaaaaa aaaaaaaaa       1859
//
'''
eb = io.StringIO(eb_str)
dna_seq = DNA.read(eb)
embl_single_record = get_data_path('embl_single_record')
dna_seq = DNA.read(embl_single_record)

# class SnifferTests(TestCase):
#     def setUp(self):
#         self.positive_fps = list(map(get_data_path, [
#             'genbank_5_blanks_start_of_file',
#             'genbank_single_record_upper',
#             'genbank_single_record_lower',
#             'genbank_multi_records']))

#         self.negative_fps = list(map(get_data_path, [
#             'empty',
#             'whitespace_only',
#             'genbank_6_blanks_start_of_file',
#             'genbank_w_beginning_whitespace',
#             'genbank_missing_locus_name']))

#     def test_positives(self):
#         for fp in self.positive_fps:
#             self.assertEqual(_genbank_sniffer(fp), (True, {}))

#     def test_negatives(self):
#         for fp in self.negative_fps:
#             self.assertEqual(_genbank_sniffer(fp), (False, {}))


# class GenBankIOTests(TestCase):
#     # parent class to set up test data for the child class
#     def setUp(self):
#         # test locus line
#         self.locus = (
#             (['LOCUS       NC_005816   9609 bp   '
#               'DNA   circular   CON   07-FEB-2015'],
#              {'division': 'CON', 'mol_type': 'DNA', 'shape': 'circular',
#               'locus_name': 'NC_005816', 'date': '07-FEB-2015',
#               'unit': 'bp', 'size': 9609}),
#             (['LOCUS       SCU49845   5028 bp   '
#               'DNA      PLN   21-JUN-1999'],
#              {'division': 'PLN', 'mol_type': 'DNA', 'shape': None,
#              'locus_name': 'SCU49845', 'date': '21-JUN-1999',
#               'unit': 'bp', 'size': 5028}),
#             (['LOCUS       NP_001832   360 aa      '
#               'linear   PRI   18-DEC-2001'],
#              {'division': 'PRI', 'mol_type': None, 'shape': 'linear',
#               'locus_name': 'NP_001832', 'date': '18-DEC-2001',
#               'unit': 'aa', 'size': 360}))

#         # test single record and read uppercase sequence
#         self.single_upper_fp = get_data_path('genbank_single_record_upper')
#         self.single_lower_fp = get_data_path('genbank_single_record_lower')
#         self.single = (
#             'GSREILDFK',
#             {'LOCUS': {'date': '23-SEP-1994',
#                        'division': 'BCT',
#                        'locus_name': 'AAB29917',
#                        'mol_type': None,
#                        'shape': 'linear',
#                        'size': 9,
#                        'unit': 'aa'}},
#             None,
#             Protein)

#         self.single_rna_fp = get_data_path('genbank_single_record')
#         imd = IntervalMetadata(63)
#         imd.add([(0, 63)],
#                 [(False, False)],
#                 {'db_xref': '"taxon:562"',
#                  'mol_type': '"mRNA"',
#                  'organism': '"Escherichia coli"',
#                  'type': 'source',
#                  'strand': '+',
#                  '__location': '1..63'})
#         imd.add([(0, 63)],
#                 [(False, True)],
#                 {'phase': 0,
#                  'db_xref': ['"GI:145230"', '"taxon:562"', '"taxon:561"'],
#                  '__location': '1..>63',
#                  'strand': '+',
#                  'note': '"alkaline phosphatase signal peptide"',
#                  'protein_id': '"AAA23431.1"',
#                  'transl_table': '11',
#                  'translation': '"MKQSTIALAVLPLLFTPVTKA"',
#                  'type': 'CDS'})
#         self.single_rna = (
#             'gugaaacaaagcacuauugcacuggcugucuuaccguuacuguuuaccccugugacaaaagcc',
#             {'ACCESSION': 'M14399',
#              'COMMENT': 'Original source text: E.coli, cDNA to mRNA.',
#              'DEFINITION': "alkaline phosphatase signal mRNA, 5' end.",
#              'KEYWORDS': 'alkaline phosphatase; signal peptide.',
#              'LOCUS': {'date': '26-APR-1993',
#                        'division': 'BCT',
#                        'locus_name': 'ECOALKP',
#                        'mol_type': 'mRNA',
#                        'shape': 'linear',
#                        'size': 63,
#                        'unit': 'bp'},
#              'SOURCE': {'ORGANISM': 'Escherichia coli',
#                         'taxonomy': 'Bacteria; Proteobacteria; '
#                         'Gammaproteobacteria; Enterobacteriales; '
#                         'Enterobacteriaceae; Escherichia.'},
#              'VERSION': 'M14399.1  GI:145229'},
#             imd,
#             RNA)

#         # test:
#         # 1. multiple records in one file
#         # 2. lowercase sequence
#         # 3. DNA, RNA, Protein type
#         # 4. variation of formats
#         self.multi_fp = get_data_path('genbank_multi_records')
#         imd_pro = IntervalMetadata(9)
#         imd_pro.add([(0, 9)], [(False, False)],
#                     {'organism': '"Bacteria"',
#                      'type': 'source',
#                      'strand': '+',
#                      '__location': '1..9'},)
#         imd_pro.add([(0, 9)], [(False, True)],
#                     {'__location': '1..>9',
#                      'product': '"L-carnitine amidase"',
#                      'strand': '+',
#                      'type': 'Protein'})
#         imd_dna = IntervalMetadata(9)
#         imd_dna.add([(0, 9)], [(False, False)],
#                     {'country': '"Brazil: Parana, Paranavai"',
#                      'type': 'source',
#                      'strand': '+',
#                      '__location': '1..9',
#                      'environmental_sample': ''})
#         imd_dna.add([(1, 8)], [(True, True)],
#                     {'__location': 'complement(<2..>8)',
#                      'product': '"16S ribosomal RNA"',
#                      'strand': '-',
#                      'type': 'rRNA'})

#         self.multi = (
#             ('gsreildfk',
#              {'ACCESSION': 'AAB29917',
#               'COMMENT': 'Method: direct peptide sequencing.',
#               'DBSOURCE': 'accession AAB29917.1',
#               'DEFINITION': 'L-carnitine amidase {N-terminal}',
#               'KEYWORDS': '.',
#               'LOCUS': {'date': '23-SEP-1994',
#                         'division': 'BCT',
#                         'locus_name': 'AAB29917',
#                         'mol_type': None,
#                         'shape': 'linear',
#                         'size': 9,
#                         'unit': 'aa'},
#               'REFERENCE': [{'AUTHORS': 'Joeres,U. and Kula,M.R.',
#                              'JOURNAL': 'AMB 40 (5), 606-610 (1994)',
#                              'PUBMED': '7764422',
#                              'REFERENCE': '1  (residues 1 to 9)',
#                              'REMARK': 'from the original journal article.',
#                              'TITLE': 'a microbial L-carnitine amidase'},
#                             {'AUTHORS': 'Joeres,U. and Kula,M.R.',
#                              'JOURNAL': 'AMB 40 (5), 606-610 (1994)',
#                              'PUBMED': '7764422',
#                              'REFERENCE': '1  (residues 1 to 9)',
#                              'TITLE': 'a microbial L-carnitine amidase'}],
#               'SOURCE': {'ORGANISM': 'Bacteria',
#                          'taxonomy': 'Unclassified.'},
#               'VERSION': 'AAB29917.1  GI:545426'},
#              imd_pro,
#              Protein),

#             ('catgcaggc',
#              {'ACCESSION': 'HQ018078',
#               'DEFINITION': 'Uncultured Xylanimonas sp.16S, partial',
#               'KEYWORDS': 'ENV.',
#               'LOCUS': {'date': '29-AUG-2010',
#                         'division': 'ENV',
#                         'locus_name': 'HQ018078',
#                         'mol_type': 'DNA',
#                         'shape': 'linear',
#                         'size': 9,
#                         'unit': 'bp'},
#               'SOURCE': {'ORGANISM': 'uncultured Xylanimonas sp.',
#                          'taxonomy': 'Bacteria; Actinobacteria; '
#                          'Micrococcales; Promicromonosporaceae; '
#                          'Xylanimonas; environmental samples.'},
#               'VERSION': 'HQ018078.1  GI:304421728'},
#              imd_dna,
#              DNA))


# class ReaderTests(GenBankIOTests):
#     def test_parse_reference(self):
#         lines = '''
# REFERENCE   1  (bases 1 to 154478)
#   AUTHORS   Sato,S., Nakamura,Y., Kaneko,T., and Tabata,S.
#   TITLE     Complete structure of the chloroplast genome of
#             Arabidopsis thaliana
#   JOURNAL   DNA Res. 6 (5), 283-290 (1999)
#    PUBMED   10574454'''.split('\n')

#         exp = {'AUTHORS': 'Sato,S., Nakamura,Y., Kaneko,T., and Tabata,S.',
#                'JOURNAL': 'DNA Res. 6 (5), 283-290 (1999)',
#                'PUBMED': '10574454',
#                'REFERENCE': '1  (bases 1 to 154478)',
#                'TITLE': ('Complete structure of the chloroplast genome of'
#                          ' Arabidopsis thaliana')}
#         self.assertEqual(_parse_reference(lines), exp)

#     def test_parse_locus(self):
#         for serialized, parsed in self.locus:
#             self.assertEqual(_parse_locus(serialized), parsed)

#     def test_parse_locus_invalid(self):
#         lines = [
#             # missing unit
#             ['LOCUS       NC_005816               9609 '
#              '    DNA     circular CON 07-FEB-2015'],
#             # missing division
#             ['LOCUS       SCU49845     5028 bp'
#              '    DNA                    21-JUN-1999'],
#             # wrong date format
#             ['LOCUS       NP_001832                360 aa'
#              '            linear   PRI 2001-12-18']]
#         for line in lines:
#             with self.assertRaisesRegex(GenBankFormatError,
#                                         'Could not parse the LOCUS line:.*'):
#                 _parse_locus(line)

#     def test_genbank_to_generator_single(self):
#         # test single record and uppercase sequence
#         for c in [Sequence, Protein]:
#             obs = next(_genbank_to_generator(
#                 self.single_upper_fp, constructor=c))
#             exp = c(self.single[0], metadata=self.single[1],
#                     positional_metadata=self.single[2])
#             self.assertEqual(exp, obs)

#     def test_genbank_to_generator(self):
#         for i, obs in enumerate(_genbank_to_generator(self.multi_fp)):
#             seq, md, imd, constructor = self.multi[i]
#             exp = constructor(seq, metadata=md, lowercase=True,
#                               interval_metadata=imd)
#             self.assertEqual(exp, obs)

#     def test_genbank_to_sequence(self):
#         for i, exp in enumerate(self.multi):
#             obs = _genbank_to_sequence(self.multi_fp, seq_num=i+1)
#             exp = Sequence(exp[0], metadata=exp[1], lowercase=True,
#                            interval_metadata=exp[2])
#             self.assertEqual(exp, obs)

#     def test_genbank_to_rna(self):
#         seq, md, imd, constructor = self.single_rna
#         obs = _genbank_to_rna(self.single_rna_fp)
#         exp = constructor(seq, metadata=md,
#                           lowercase=True, interval_metadata=imd)

#         self.assertEqual(exp, obs)

#     def test_genbank_to_dna(self):
#         i = 1
#         exp = self.multi[i]
#         obs = _genbank_to_dna(self.multi_fp, seq_num=i+1)
#         exp = DNA(exp[0], metadata=exp[1], lowercase=True,
#                   interval_metadata=exp[2])

#         self.assertEqual(exp, obs)

#     def test_genbank_to_protein(self):
#         i = 0
#         exp = self.multi[i]
#         obs = _genbank_to_protein(self.multi_fp, seq_num=i+1)
#         exp = Protein(exp[0], metadata=exp[1],
#                       lowercase=True, interval_metadata=exp[2])
#         self.assertEqual(exp, obs)


# class WriterTests(GenBankIOTests):
#     def test_serialize_locus(self):
#         for serialized, parsed in self.locus:
#             self.assertEqual(
#                 _serialize_locus('LOCUS', parsed), serialized[0] + '\n')

#     def test_generator_to_genbank(self):
#         seq, md, imd, constructor = self.single
#         obj = constructor(seq, md, interval_metadata=imd)
#         with io.StringIO() as fh:
#             _generator_to_genbank([obj], fh)
#             obs = fh.getvalue()

#         with open(self.single_lower_fp) as fh:
#             exp = fh.read()

#         self.assertEqual(obs, exp)

#     def test_sequence_to_genbank(self):
#         with io.StringIO() as fh:
#             for i, (seq, md, imd, constructor) in enumerate(self.multi):
#                 obj = Sequence(seq, md, interval_metadata=imd, lowercase=True)
#                 _sequence_to_genbank(obj, fh)
#             obs = fh.getvalue()

#         with open(self.multi_fp) as fh:
#             exp = fh.read()

#         self.assertEqual(obs, exp)

#     def test_dna_protein_to_genbank(self):
#         writers = [_protein_to_genbank,
#                    _dna_to_genbank]
#         with io.StringIO() as fh:
#             for i, (seq, md, imd, constructor) in enumerate(self.multi):
#                 obj = constructor(
#                     seq, md, interval_metadata=imd, lowercase=True)
#                 writers[i](obj, fh)
#             obs = fh.getvalue()

#         with open(self.multi_fp) as fh:
#             exp = fh.read()

#         self.assertEqual(obs, exp)

#     def test_rna_to_genbank(self):
#         with io.StringIO() as fh:
#             seq, md, imd, constructor = self.single_rna
#             obj = constructor(seq, md, interval_metadata=imd, lowercase=True)
#             _rna_to_genbank(obj, fh)
#             obs = fh.getvalue()

#         with open(self.single_rna_fp) as fh:
#             exp = fh.read()

#         self.assertEqual(obs, exp)


# class RoundtripTests(GenBankIOTests):
#     def test_roundtrip_generator(self):
#         with io.StringIO() as fh:
#             _generator_to_genbank(_genbank_to_generator(self.multi_fp), fh)
#             obs = fh.getvalue()

#         with open(self.multi_fp) as fh:
#             exp = fh.read()

#         self.assertEqual(obs, exp)

#     def test_roundtrip_rna(self):
#         with io.StringIO() as fh:
#             _rna_to_genbank(_genbank_to_rna(self.single_rna_fp), fh)
#             obs = fh.getvalue()

#         with open(self.single_rna_fp) as fh:
#             exp = fh.read()

#         self.assertEqual(obs, exp)

#     def test_roundtrip_dna(self):
#         with io.StringIO() as fh:
#             _dna_to_genbank(_genbank_to_dna(self.single_rna_fp), fh)
#             obs = fh.getvalue()

#         with open(self.single_rna_fp) as fh:
#             exp = fh.read()

#         self.assertEqual(obs, exp)

#     def test_roundtrip_protein(self):
#         with io.StringIO() as fh:
#             _protein_to_genbank(_genbank_to_protein(self.single_lower_fp), fh)
#             obs = fh.getvalue()

#         with open(self.single_lower_fp) as fh:
#             exp = fh.read()

#         self.assertEqual(obs, exp)

#     def test_roundtrip_sequence(self):
#         with io.StringIO() as fh:
#             _sequence_to_genbank(_genbank_to_sequence(self.single_rna_fp), fh)
#             obs = fh.getvalue()

#         with open(self.single_rna_fp) as fh:
#             exp = fh.read()

#         self.assertEqual(obs, exp)


# if __name__ == '__main__':
#     main()
