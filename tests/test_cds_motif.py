from os import path
import unittest
from Bio import SeqIO
from secmet.record import Record

filename = 'nisin.gbk'
filetype = 'genbank'

class TestCDS_motifFeature(unittest.TestCase):
    def get_testfile(self):
        """File path for testing"""
        return path.join(path.dirname(__file__), 'data', filename)

    def test_CDS_motifFeature_members(self):
        """Check if all the qualifiers are properly stored in CDS_motifFeature"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_cds_motifs = [i for i in bp_rec.features if i.type == 'CDS_motif']
        mod_cds_motifs = rec.get_CDS_motifs()
        for bp_motif, mod_motif in zip(bp_cds_motifs, mod_cds_motifs):
            for key, value in bp_motif.qualifiers.items():
                if key == 'note':
                    #note is modified to notes in secmet
                    self.assertEqual(bp_motif.qualifiers['note'], mod_motif.notes)
                else:
                    if value is not None and value:
                        if not hasattr(mod_motif, key):
                            raise AttributeError("%s is not a member of CDS_motifFeature"%key)
                        self.assertEqual(str(value[0]), str(getattr(mod_motif, key)))
