from os import path
import unittest
from Bio import SeqIO
from secmet.record import Record

filename = 'nisin.gbk'
filetype = 'genbank'

class TestDomains(unittest.TestCase):
    def get_testfile(self):
        """File path for testing"""
        return path.join(path.dirname(__file__), 'data', filename)

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
                    if key != 'note':
                        self.assertEqual(value, mod_generic.get_qualifier(key))
                    else:
                        #note is modified to notes in secmet
                        self.assertEqual(bp_generic.qualifiers['note'], mod_generic.notes)
