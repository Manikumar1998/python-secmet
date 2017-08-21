from os import path
import unittest
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from secmet.record import Record, CDS_motifFeature

filename = 'Y16952.3.final.gbk'
filetype = 'genbank'

class TestCDS_motifFeature(unittest.TestCase):
    def get_testfile(self):
        """File path for testing"""
        return path.join(path.dirname(__file__), 'data', filename)

    def BioFeature(self):
        biofeature = SeqFeature(location=FeatureLocation(10, 100))
        biofeature.qualifiers = {'locus_tag': ['fake_locus_tag'], 'translation': ['fake_translation'], 'aSTool': ['fake_aSTool'],\
                                 'motif': ['fake_motif'], 'asDomain_id': ['fake_asDomain_id'], 'detection': ['fake_detection'],\
                                 'database': ['fake_database'], 'label': ['fake_label'], 'unknown_qualifier': ['fake_qualifier'],\
                                 'note': ['fake_notes'], 'aSProdPred': ['fake_aSProdPred'], 'aSASF_choice':['fake_aSAF_choice'],\
                                 'aSASF_note': ['fake_aSASF_note'], 'aSASF_choice': ['fake_aSASF_choice'], 'aSASF_scaffold': ['fake_aSASF_scaffold'],\
                                 'aSASF_prediction': ['fake_aSASF_prediction']}
        return biofeature

    def test_CDS_motifFeature_members(self):
        """Check if all the qualifiers are properly stored in CDS_motifFeature"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_cds_motifs = [i for i in bp_rec.features if i.type == 'CDS_motif']
        mod_cds_motifs = rec.get_CDS_motifs()
        #Segregate out qualifiers that are stored in list form
        qualifiers_as_list = ['note', 'aSASF_choice', 'aSASF_note', 'aSASF_scaffold', \
                              'aSASF_prediction', 'aSProdPred']
        for bp_motif, mod_motif in zip(bp_cds_motifs, mod_cds_motifs):
            for key, value in bp_motif.qualifiers.items():
                if value:
                    if key not in qualifiers_as_list:
                        if not hasattr(mod_motif, key):
                            raise AttributeError("%s is not a member of CDS_motifFeature"%key)
                        #score and evalue are numbers
                        if key in ['score', 'evalue']:
                            self.assertEqual(float(value[0]), float(getattr(mod_motif, key)))
                        else:
                            self.assertEqual(str(value[0]), str(getattr(mod_motif, key)))
                    else:
                        if key == 'note':
                            #note is modified to notes in secmet
                            self.assertEqual(value, mod_motif.notes)
                        else:
                            self.assertEqual(value, getattr(mod_motif, key))
        cdsmotif = CDS_motifFeature(FeatureLocation(1, 10))
        #score, evalue should be numbers
        try:
            cdsmotif.score = '-a50'
        except ValueError:
            pass
        try:
            cdsmotif.evalue = 'a5.50E-08'
        except ValueError:
            pass

        #If valid qualifiers and values are added, We shouldn't get an error
        try:
            cdsmotif.score = '-50'
            cdsmotif.evalue = '5.50E-08'
        except:
            raise RuntimeError('Secmet unable to add valid qualifiers')

    def test_BioFeature_to_CDS_motifFeauture(self):
        biofeature = self.BioFeature()
        cdsmotif_feature = CDS_motifFeature(feature=biofeature)
        self.assertEqual(str(cdsmotif_feature.location), str(FeatureLocation(10, 100)))
        self.assertEqual(cdsmotif_feature.type, 'CDS_motif')
        self.assertEqual(cdsmotif_feature.locus_tag, 'fake_locus_tag')
        self.assertEqual(cdsmotif_feature.translation, 'fake_translation')
        self.assertEqual(cdsmotif_feature.label, 'fake_label')
        self.assertEqual(cdsmotif_feature.aSTool, 'fake_aSTool')
        self.assertEqual(cdsmotif_feature.detection, 'fake_detection')
        self.assertEqual(cdsmotif_feature.database, 'fake_database')
        self.assertEqual(cdsmotif_feature.asDomain_id, 'fake_asDomain_id')
        self.assertEqual(cdsmotif_feature.motif, 'fake_motif')
        self.assertEqual(cdsmotif_feature.notes, ['fake_notes'])
        self.assertEqual(cdsmotif_feature.aSProdPred, ['fake_aSProdPred'])
        self.assertEqual(cdsmotif_feature.aSASF_note, ['fake_aSASF_note'])
        self.assertEqual(cdsmotif_feature.aSASF_scaffold, ['fake_aSASF_scaffold'])
        self.assertEqual(cdsmotif_feature.aSASF_choice, ['fake_aSASF_choice'])
        self.assertEqual(cdsmotif_feature.aSASF_prediction, ['fake_aSASF_prediction'])
        self.assertEqual(cdsmotif_feature._qualifiers['unknown_qualifier'], ['fake_qualifier'])
        self.assertEqual(repr(cdsmotif_feature), repr(cdsmotif_feature.to_biopython()[0]))
        self.assertIsInstance(cdsmotif_feature.to_biopython()[0], SeqFeature)
