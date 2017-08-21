from os import path
import unittest
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from secmet.record import Record, aSDomain, PFAM_domain

filename = 'Y16952.3.final.gbk'
filetype = 'genbank'

class TestDomains(unittest.TestCase):
    def get_testfile(self):
        """File path for testing"""
        return path.join(path.dirname(__file__), 'data', filename)

    def BioFeature(self):
        biofeature = SeqFeature(location=FeatureLocation(10, 100))
        biofeature.qualifiers = {'locus_tag': ['fake_locus_tag'], 'translation': ['fake_translation'], 'aSTool': ['fake_aSTool'], \
                                 'domain': ['fake_domain'], 'asDomain_id': ['fake_asDomain_id'], 'detection': ['fake_detection'], \
                                 'database': ['fake_database'], 'label': ['fake_label'], 'unknown_qualifier': ['fake_qualifier'], \
                                 'db_xref': ['fake_db_xref'], 'note': ['fake_notes'], 'aSProdPred': ['fake_aSProdPred'], \
                                 'domain_subtype': ['fake_domain_subtype'], 'aSASF_choice':['fake_aSAF_choice'], 'aSASF_note': ['fake_aSASF_note'], \
                                 'aSASF_choice': ['fake_aSASF_choice'], 'aSASF_scaffold': ['fake_aSASF_scaffold'], 'aSASF_prediction': ['fake_aSASF_prediction'], \
                                 'specificity': ['fake_specificity']}
        return biofeature

    def test_aSDomain(self):
        """Check if all the qualifiers are properly stored in aSDomain"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_asdomains = [i for i in bp_rec.features if i.type == 'aSDomain']
        mod_asdomains = rec.get_aSDomains()
        #Segregate out qualifiers that are stored in list form
        qualifiers_as_list = ['note', 'specificity', 'aSASF_choice', 'aSASF_note', \
                              'aSASF_scaffold', 'aSASF_prediction', 'aSProdPred']
        for bp_asdomain, mod_asdomain in zip(bp_asdomains, mod_asdomains):
            for key, value in bp_asdomain.qualifiers.items():
                if value:
                    if key not in qualifiers_as_list:
                        if not hasattr(mod_asdomain, key):
                            raise AttributeError('%s is not a member of aSDomain'%key)
                        if key in ['score', 'evalue']:
                            self.assertEqual(float(value[0]), float(getattr(mod_asdomain, key)))
                        else:
                            self.assertEqual(str(value[0]), str(getattr(mod_asdomain, key)))
                    else:
                        if key == 'note':
                            #note is modified to notes in secmet
                            self.assertEqual(value, mod_asdomain.notes)
                        else:
                            self.assertEqual(value, getattr(mod_asdomain, key))

        asdomain = aSDomain(FeatureLocation(1, 10))
        #score, evalue should be numbers
        try:
            asdomain.score = '-a50'
        except ValueError:
            pass
        try:
            asdomain.evalue = 'a5.50E-08'
        except ValueError:
            pass

        #If valid qualifiers and values are added, We shouldn't get an error
        try:
            asdomain.score = '-50'
            asdomain.evalue = '5.50E-08'
        except:
            raise RuntimeError('Secmet unable to add valid qualifiers')

    def test_BioFeature_to_aSDomain(self):
        biofeature = self.BioFeature()
        asdomain_feature = aSDomain(feature=biofeature)
        self.assertEqual(str(asdomain_feature.location), str(FeatureLocation(10, 100)))
        self.assertEqual(asdomain_feature.type, 'aSDomain')
        self.assertEqual(asdomain_feature.locus_tag, 'fake_locus_tag')
        self.assertEqual(asdomain_feature.translation, 'fake_translation')
        self.assertEqual(asdomain_feature.label, 'fake_label')
        self.assertEqual(asdomain_feature.detection, 'fake_detection')
        self.assertEqual(asdomain_feature.database, 'fake_database')
        self.assertEqual(asdomain_feature.asDomain_id, 'fake_asDomain_id')
        self.assertEqual(asdomain_feature.domain, 'fake_domain')
        self.assertEqual(asdomain_feature.domain_subtype, 'fake_domain_subtype')
        self.assertEqual(asdomain_feature.specificity, ['fake_specificity'])
        self.assertEqual(asdomain_feature.notes, ['fake_notes'])
        self.assertEqual(asdomain_feature.aSProdPred, ['fake_aSProdPred'])
        self.assertEqual(asdomain_feature.aSASF_note, ['fake_aSASF_note'])
        self.assertEqual(asdomain_feature.aSASF_scaffold, ['fake_aSASF_scaffold'])
        self.assertEqual(asdomain_feature.aSASF_choice, ['fake_aSASF_choice'])
        self.assertEqual(asdomain_feature.aSASF_prediction, ['fake_aSASF_prediction'])
        self.assertEqual(asdomain_feature._qualifiers['unknown_qualifier'], ['fake_qualifier'])
        self.assertEqual(repr(asdomain_feature), repr(asdomain_feature.to_biopython()[0]))
        self.assertIsInstance(asdomain_feature.to_biopython()[0], SeqFeature)

    def test_PFAM_domain(self):
        """Check if all the qualifiers are properly stored in PFAM_domain"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_pfams = [i for i in bp_rec.features if i.type == 'PFAM_domain']
        mod_pfams = rec.get_PFAM_domains()
        #Segregate out qualifiers that are stored in list form
        qualifiers_as_list = ['note', 'db_xref', 'aSASF_choice', 'aSASF_note', \
                              'aSASF_scaffold', 'aSASF_prediction', 'aSProdPred']
        for bp_pfam, mod_pfam in zip(bp_pfams, mod_pfams):
            for key, value in bp_pfam.qualifiers.items():
                if value:
                    if key not in qualifiers_as_list:
                        if not hasattr(mod_pfam, key):
                            raise AttributeError('%s is not a member of PFAM_domain'%key)
                        if key in ['score', 'evalue']:
                            self.assertEqual(float(value[0]), float(getattr(mod_pfam, key)))
                        else:
                            self.assertEqual(str(value[0]), str(getattr(mod_pfam, key)))
                    else:
                        if key == 'note':
                            #note is modified to notes in secmet
                            self.assertEqual(value, mod_pfam.notes)
                        else:
                            self.assertEqual(value, getattr(mod_pfam, key))

        pfam = PFAM_domain(FeatureLocation(1, 10))
        #score, evalue should be numbers
        try:
            pfam.score = '-a50'
        except ValueError:
            pass
        try:
            pfam.evalue = 'a5.50E-08'
        except ValueError:
            pass

        #If valid qualifiers and values are added, We shouldn't get an error
        try:
            pfam.score = '-50'
            pfam.evalue = '5.50E-08'
        except:
            raise RuntimeError('Secmet unable to add valid qualifiers')


    def test_BioFeature_to_PFAM_domain(self):
        biofeature = self.BioFeature()
        pfam_feature = PFAM_domain(feature=biofeature)
        self.assertEqual(str(pfam_feature.location), str(FeatureLocation(10, 100)))
        self.assertEqual(pfam_feature.type, 'PFAM_domain')
        self.assertEqual(pfam_feature.locus_tag, 'fake_locus_tag')
        self.assertEqual(pfam_feature.translation, 'fake_translation')
        self.assertEqual(pfam_feature.label, 'fake_label')
        self.assertEqual(pfam_feature.aSTool, 'fake_aSTool')
        self.assertEqual(pfam_feature.detection, 'fake_detection')
        self.assertEqual(pfam_feature.database, 'fake_database')
        self.assertEqual(pfam_feature.asDomain_id, 'fake_asDomain_id')
        self.assertEqual(pfam_feature.domain, 'fake_domain')
        self.assertEqual(pfam_feature.db_xref, ['fake_db_xref'])
        self.assertEqual(pfam_feature.notes, ['fake_notes'])
        self.assertEqual(pfam_feature.aSProdPred, ['fake_aSProdPred'])
        self.assertEqual(pfam_feature.aSASF_note, ['fake_aSASF_note'])
        self.assertEqual(pfam_feature.aSASF_scaffold, ['fake_aSASF_scaffold'])
        self.assertEqual(pfam_feature.aSASF_choice, ['fake_aSASF_choice'])
        self.assertEqual(pfam_feature.aSASF_prediction, ['fake_aSASF_prediction'])
        self.assertEqual(pfam_feature._qualifiers['unknown_qualifier'], ['fake_qualifier'])
        self.assertEqual(repr(pfam_feature), repr(pfam_feature.to_biopython()[0]))
        self.assertIsInstance(pfam_feature.to_biopython()[0], SeqFeature)
