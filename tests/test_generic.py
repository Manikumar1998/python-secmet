from os import path
import unittest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from secmet.record import Record, GenericFeature

filename = 'Y16952.3.final.gbk'
filetype = 'genbank'

class TestDomains(unittest.TestCase):

    def get_testfile(self):
        """File path for testing"""
        return path.join(path.dirname(__file__), 'data', filename)

    def BioFeature(self):
        biofeature = SeqFeature(location=FeatureLocation(10, 100))
        biofeature.type = 'FAKE_BIO_FEATURE'
        biofeature.qualifiers = {'locus_tag': ['fake_locus_tag'], 'translation': ['fake_translation'], \
                                 'gene': ['fake_gene'], 'name': ['fake_name'], 'seq': [Seq('FAKE')], \
                                 'description': ['fake_description'], 'sec_met': ['fake_sec_met'], \
                                 'note': ['fake_notes'], 'unknown_qualifier': ['fake_qualifier']}
        return biofeature

    def test_GenericFeature(self):
        """Check if all the qualifiers are properly stored in GenericFeature"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        non_generic_features = ['CDS', 'cluster', 'CDS_motif', 'aSDomain', 'PFAM_domain']
        bp_generics = [i for i in bp_rec.features if i.type not in non_generic_features]
        mod_generics = rec.get_generics()
        self.assertEqual(len(bp_generics), len(mod_generics))
        for bp_generic, mod_generic in zip(bp_generics, mod_generics):
            for key, value in bp_generic.qualifiers.items():
                if value is not None and value:
                    if key == 'note':
                        #note is modified to notes in secmet
                        self.assertEqual(bp_generic.qualifiers['note'], mod_generic.notes)
                    else:
                        self.assertEqual(value, mod_generic.get_qualifier(key))

    def test_add_qualifier(self):
        """Test adding a new GenericFeature"""
        #GenericFeature should be initialised with valid location and type
        try:
            new_generic = GenericFeature()
        except TypeError:
            pass
        try:
            new_generic = GenericFeature(FeatureLocation(10, 100), 20) #Invalid type
        except ValueError:
            pass
        new_generic = GenericFeature(FeatureLocation(1, 100), 'FAKE')

        #qualifiers should be strings and their values should be either strings or list of strings
        try:
            new_generic.add_qualifier(10, (20, 30))
        except TypeError:
            pass

        #If the formats are valid shouldn't get any error
        try:
            new_generic.add_qualifier('Some string1', 'Fake_value1')
            new_generic.add_qualifier('Some string2', ['Fake_value2'])
        except:
            raise RuntimeError('Secmet unable to add valid qualifiers')

        #score, evalue and probability should be numbers
        try:
            new_generic.add_qualifier('score', '-a50')
        except ValueError:
            pass
        try:
            new_generic.add_qualifier('evalue', 'a5.50E-08')
        except ValueError:
            pass
        try:
            new_generic.add_qualifier('probability', 'a0.5')
        except ValueError:
            pass

        #If valid qualifiers and values are added, We shouldn't get an error
        try:
            new_generic.add_qualifier('score', '-50')
            new_generic.add_qualifier('evalue', '5.50E-08')
            new_generic.add_qualifier('probability', '0.5')
        except:
            raise RuntimeError('Secmet unable to add valid qualifiers')

        #If GenericFeature has the qualifier as member, the member should get initialised
        new_generic.add_qualifier('locus_tag', 'FAKE_TAG')
        self.assertEqual(new_generic.locus_tag, 'FAKE_TAG')

        #If GenericFeature has the qualifier as memeber and if the memeber is a list, the new qualifier value should get appended
        new_generic.add_qualifier('sec_met', 'FAKE_sec_met1')
        self.assertEqual(new_generic.sec_met, ['FAKE_sec_met1'])
        new_generic.add_qualifier('sec_met', 'FAKE_sec_met2')
        self.assertEqual(new_generic.sec_met, ['FAKE_sec_met1', 'FAKE_sec_met2'])
        #If the qualifier is a list and value is also a list, qualifier should extend to new values
        new_generic.add_qualifier('sec_met', ['FAKE_sec_met3'])
        self.assertEqual(new_generic.sec_met, ['FAKE_sec_met1', 'FAKE_sec_met2', 'FAKE_sec_met3'])

        #If GenericFeature doesn't contain the qualifier as its member, the member gets stored in _qualifiers
        new_generic.add_qualifier('fake_qualifier1', 'FAKE1')
        self.assertEqual(new_generic._qualifiers['fake_qualifier1'], ['FAKE1'])

        #If a new value is added for existing qualifier it should get appended
        new_generic.add_qualifier('fake_qualifier1', 'FAKE2')
        self.assertEqual(new_generic._qualifiers['fake_qualifier1'], ['FAKE1', 'FAKE2'])

        #If the qualifier value is a list, then it should store it as a list only
        new_generic.add_qualifier('fake_qualifier2', ['FAKE1'])

    def test_convert_Biofeature_to_GenericFeature(self):
        """Test the convesion of BioFeature to GenericFeature"""
        biofeature = self.BioFeature()
        generic_feature = GenericFeature(feature=biofeature)
        self.assertEqual(str(generic_feature.location), str(FeatureLocation(10, 100)))
        self.assertEqual(generic_feature.type, 'FAKE_BIO_FEATURE')
        self.assertEqual(generic_feature.locus_tag, 'fake_locus_tag')
        self.assertEqual(generic_feature.translation, 'fake_translation')
        self.assertEqual(generic_feature.gene, 'fake_gene')
        self.assertEqual(generic_feature.name, 'fake_name')
        self.assertEqual(generic_feature.seq, Seq('FAKE'))
        self.assertEqual(generic_feature.description, 'fake_description')
        self.assertEqual(generic_feature.sec_met, ['fake_sec_met'])
        self.assertEqual(generic_feature.notes, ['fake_notes'])
        self.assertEqual(generic_feature._qualifiers['unknown_qualifier'], ['fake_qualifier'])
        self.assertEqual(repr(generic_feature), repr(generic_feature.to_biopython()[0]))
        self.assertIsInstance(generic_feature.to_biopython()[0], SeqFeature)

    def test_get_qualifier(self):
        """Test the get_qualifier method of GenericFeature"""
        biofeature = self.BioFeature()
        generic_feature = GenericFeature(feature=biofeature)
        #getting qualifiers should be case insensitive
        self.assertEqual(generic_feature.get_qualifier('LoCUs_tAg'), ['fake_locus_tag'])
        #Upper case
        generic_feature.add_qualifier('SOME_UPPER_CASE_QUALIFIER', 'FAKE')
        self.assertEqual(generic_feature.get_qualifier('some_upper_case_qualifier'), ['FAKE'])
        #Lower case
        generic_feature.add_qualifier('some_lower_case_qualifier', 'FAKE')
        self.assertEqual(generic_feature.get_qualifier('SOME_LOWER_CASE_QUALIFIER'), ['FAKE'])

        #If qualifier is not present, [] should be returned
        self.assertEqual(generic_feature.get_qualifier('Absent_qualifier'), [])
