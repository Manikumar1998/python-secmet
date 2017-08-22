from os import path
import unittest
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from secmet.record import Record, CDSFeature, SecMetQualifier, SecMetResult

filename = 'Y16952.3.final.gbk'
filetype = 'genbank'

class FakeResult(object):
    """A FakeResult to test SecMetResult"""
    def __init__(self):
        """Initialise members with fake values"""
        self.query_id = 'fake_id'
        self.evalue = '10000'
        self.bitscore = '10000'

class TestCDSFeature(unittest.TestCase):
    def get_testfile(self):
        """File path for testing"""
        return path.join(path.dirname(__file__), 'data', filename)

    def BioFeature(self):
        biofeature = SeqFeature(location=FeatureLocation(10, 100))
        biofeature.qualifiers = {'locus_tag': ['fake_locus_tag'], 'translation': ['fake_translation'], \
                                 'gene': ['fake_gene'], 'product': ['fake_product'], 'protein_id': ['fake_protein_id'], \
                                 'transl_table': ['fake_transl_table'], 'source': ['fake_source'], 'db_xref': ['fake_db_xref'],\
                                 'EC_number': ['fake_EC_number'], 'note': ['fake_notes'], 'aSProdPred': ['fake_aSProdPred'], \
                                 'aSASF_choice':['fake_aSAF_choice'], 'aSASF_note': ['fake_aSASF_note'], 'aSASF_choice': ['fake_aSASF_choice'], \
                                 'aSASF_scaffold': ['fake_aSASF_scaffold'], 'aSASF_prediction': ['fake_aSASF_prediction'], \
                                 'sec_met_predictions': ['fake_sec_met_predictions'], 'unknown_qualifier': ['fake_qualifier']}
        return biofeature

    def test_CDSFeature_members(self):
        """Test the members of CDSFeature"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_cdss = [i for i in bp_rec.features if i.type == 'CDS']
        mod_cdss = rec.get_CDSs()
        self.assertEqual(len(bp_cdss), len(mod_cdss))
        #Segregate out qualifiers that are stored in list form
        qualifiers_as_list = ['note', 'aSProdPred', 'aSASF_choice', 'aSASF_prediction', 'db_xref', \
                              'aSASF_note', 'aSASF_scaffold', 'sec_met_predictions', 'EC_number']
        for bp_cds, mod_cds in zip(bp_cdss, mod_cdss):
            for key, value in bp_cds.qualifiers.items():
                if value:
                    if key not in qualifiers_as_list:
                        if key != 'sec_met':    #antiSMASH anyways erases all sec_met qualifiers
                            if not hasattr(mod_cds, key):
                                if not key in mod_cds._qualifiers:
                                    raise AttributeError('%s is not a member of CDSFeature'%key)
                                else:
                                    self.assertEqual(value, mod_cds._qualifiers[key])
                            else:
                                self.assertEqual(str(value[0]), str(getattr(mod_cds, key)))
                        else:
                            self.assertEqual(value, mod_cds.sec_met.as_list())
                    else:
                        if key == 'note':
                            #note is modified to notes in secmet
                            self.assertEqual(value, mod_cds.notes)
                        else:
                            self.assertEqual(value, getattr(mod_cds, key))

    def test_BioFeature_to_CDSFeature(self):
        biofeature = self.BioFeature()
        cds_feature = CDSFeature(feature=biofeature)
        self.assertEqual(str(cds_feature.location), str(FeatureLocation(10, 100)))
        self.assertEqual(cds_feature.type, 'CDS')
        self.assertEqual(cds_feature.locus_tag, 'fake_locus_tag')
        self.assertEqual(cds_feature.translation, 'fake_translation')
        self.assertEqual(cds_feature.gene, 'fake_gene')
        self.assertEqual(cds_feature.product, 'fake_product')
        self.assertEqual(cds_feature.protein_id, 'fake_protein_id')
        self.assertEqual(cds_feature.transl_table, 'fake_transl_table')
        self.assertEqual(cds_feature.source, 'fake_source')
        self.assertEqual(cds_feature.EC_number, ['fake_EC_number'])
        self.assertEqual(cds_feature.notes, ['fake_notes'])
        self.assertEqual(cds_feature.db_xref, ['fake_db_xref'])
        self.assertEqual(cds_feature.aSProdPred, ['fake_aSProdPred'])
        self.assertEqual(cds_feature.aSASF_note, ['fake_aSASF_note'])
        self.assertEqual(cds_feature.aSASF_scaffold, ['fake_aSASF_scaffold'])
        self.assertEqual(cds_feature.aSASF_choice, ['fake_aSASF_choice'])
        self.assertEqual(cds_feature.aSASF_prediction, ['fake_aSASF_prediction'])
        self.assertEqual(cds_feature.sec_met_predictions, ['fake_sec_met_predictions'])
        self.assertEqual(cds_feature._qualifiers['unknown_qualifier'], ['fake_qualifier'])
        self.assertEqual(repr(cds_feature), repr(cds_feature.to_biopython()[0]))
        self.assertIsInstance(cds_feature.to_biopython()[0], SeqFeature)
        self.assertIsInstance(cds_feature.sec_met, SecMetQualifier)

    def test_SecMetQualifier(self):
        """Test SecMetQualifier"""
        cds = CDSFeature(FeatureLocation(1, 10))
        self.assertEqual(None, cds.sec_met.clustertype)
        self.assertEqual(None, cds.sec_met.domains)
        self.assertEqual(None, cds.sec_met.kind)
        self.assertEqual(0, len(cds.sec_met))
        self.assertEqual([], cds.sec_met)
        self.assertEqual([], cds.sec_met.nrpspks)
        self.assertEqual([], cds.sec_met.asf_predictions)
        self.assertEqual([], cds.sec_met.as_list())

        cds.sec_met.clustertype = "FAKE"
        cds.sec_met.domains = ["FAKE_DOMAIN1", "FAKE_DOMAIN2"]
        cds.sec_met.kind = "FAKE"
        cds.sec_met.nrpspks = ["FAKE_NRPS/PKS Domain: "]
        cds.sec_met.asf_predictions = ['FAKE_ASF_predictions: ']
        self.assertEqual("FAKE", cds.sec_met.clustertype)
        self.assertEqual(["FAKE_DOMAIN1", "FAKE_DOMAIN2"], cds.sec_met.domains)
        self.assertEqual("FAKE", cds.sec_met.kind)
        self.assertEqual(5, len(cds.sec_met))
        expected_sec_met = ['Type: FAKE', 'Domains detected: FAKE_DOMAIN1; FAKE_DOMAIN2', 'Kind: FAKE', \
                            'FAKE_NRPS/PKS Domain: ', 'FAKE_ASF_predictions: ']
        self.assertEqual(expected_sec_met, cds.sec_met.as_list())
        self.assertEqual(["FAKE_NRPS/PKS Domain: "], cds.sec_met.nrpspks)
        self.assertEqual(['FAKE_ASF_predictions: '], cds.sec_met.asf_predictions)
        self.assertEqual(str(cds.sec_met), repr(cds.sec_met))
        with self.assertRaises(TypeError):
            #sec_met feature should be an instance of SecMetQualifier
            cds.sec_met = []
            cds.to_biopython()

        #Test the failure cases in SecMetQualifier
        with self.assertRaises(TypeError):
            #clustertype should be a string
            SecMetQualifier(clustertype=1)
        with self.assertRaises(TypeError):
            #domains should be a list
            SecMetQualifier(domains='invalid_domains_type')
        with self.assertRaises(TypeError):
            #kind should be a str
            SecMetQualifier(kind=1)

    def test_SecMetResult(self):
        """Test the SecMetResult class"""
        empty_result = SecMetResult()
        self.assertEqual(None, empty_result.query_id)
        self.assertEqual(None, empty_result.evalue)
        self.assertEqual(None, empty_result.bitscore)
        self.assertEqual(None, empty_result.nseeds)
        result = SecMetResult(FakeResult(), "fake_seeds")
        self.assertEqual('fake_id', result.query_id)
        self.assertEqual(10000.0, result.evalue)
        self.assertEqual(10000.0, result.bitscore)
        self.assertEqual('fake_seeds', result.nseeds)

        expected = "fake_id (E-value: 10000.0, bitscore: 10000.0, seeds: fake_seeds)"
        self.assertEqual(expected, repr(result), str(result))
        #Test the failure cases in SecMetResult
        result = SecMetResult()
        with self.assertRaises(ValueError):
            #evalue should be a float
            result.evalue = 'invalid_evalue'
        with self.assertRaises(ValueError):
            #bitscore should be a float
            result.bitscore = 'invalid_bitscore'
