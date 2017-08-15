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
        for bp_cds, mod_cds in zip(bp_cdss, mod_cdss):
            for key, value in bp_cds.qualifiers.items():
                if value is not None and value:
                    #aSProdPred, aSASF_choice, aSASF_choice, aSASF_note, aSASF_prediction
                    #aSASF_scaffold and sec_met_predictions are lists
                    if key not in ['aSProdPred', 'aSASF_choice', 'aSASF_note', 'aSASF_prediction', \
                                   'aSASF_scaffold', 'sec_met_predictions']:
                        if key != 'sec_met':    #antiSMASH anyways erases all sec_met qualifiers
                            if hasattr(mod_cds, key):
                                self.assertEqual(str(value[0]), str(getattr(mod_cds, key)))
                    else:
                        self.assertEqual(value, getattr(mod_cds, key))
            if 'note' in bp_cds.qualifiers:
                if bp_cds.qualifiers['note']:
                    #note is modified to notes in secmet
                    self.assertEqual(bp_cds.qualifiers['note'], mod_cds.notes)
