from os import path
import unittest
from Bio import SeqIO
from secmet.record import Record

filename = 'nisin.gbk'
filetype = 'genbank'

class TestCDSFeature(unittest.TestCase):
    def get_testfile(self):
        """File path for testing"""
        return path.join(path.dirname(__file__), 'data', filename)

    def test_CDSFeature_members(self):
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_cdss = [i for i in bp_rec.features if i.type == 'CDS']
        mod_cdss = rec.get_CDSs()
        self.assertEqual(len(bp_cdss), len(mod_cdss))
        #aSProdPred, aSASF_choice, aSASF_choice, aSASF_note, aSASF_prediction
        #aSASF_scaffold and sec_met_predictions are lists
        qualifiers_as_list = ['note', 'aSProdPred', 'aSASF_choice', 'aSASF_prediction', \
                              'aSASF_note', 'aSASF_scaffold', 'sec_met_predictions']
        for bp_cds, mod_cds in zip(bp_cdss, mod_cdss):
            for key, value in bp_cds.qualifiers.items():
                if value is not None and value:
                    if key not in qualifiers_as_list:
                        if key != 'sec_met':    #antiSMASH anyways erases all sec_met qualifiers
                            if not hasattr(mod_cds, key):
                                if not key in mod_cds._qualifiers:
                                    raise AttributeError('%s is not a member of CDSFeature'%key)
                                else:
                                    self.assertEqual(bp_cds.qualifiers[key], mod_cds._qualifiers[key])
                            else:
                                self.assertEqual(str(value[0]), str(getattr(mod_cds, key)))
                    else:
                        if key == 'note':
                            #note is modified to notes in secmet
                            self.assertEqual(value, mod_cds.notes)
                        else:
                            self.assertEqual(value, getattr(mod_cds, key))
