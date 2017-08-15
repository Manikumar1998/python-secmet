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

    def test_aSDomain(self):
        """Check if all the qualifiers are properly stored in aSDomain"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_asdomains = [i for i in bp_rec.features if i.type == 'aSDomain']
        mod_asdomains = rec.get_aSDomains()
        qualifiers_as_list = ['note', 'label', 'specificity']
        for bp_asdomain, mod_asdomain in zip(bp_asdomains, mod_asdomains):
            for key, value in bp_asdomain.qualifiers.items():
                if value is not None and value:
                    #label and specificity are lists
                    if key not in qualifiers_as_list:
                        if not hasattr(mod_asdomain, key):
                            raise AttributeError('%s is not a member of aSDomain'%key)
                        self.assertEqual(str(value[0]), str(getattr(mod_asdomain, key)))
                    else:
                        if key == 'note':
                            #note is modified to notes in secmet
                            self.assertEqual(value, mod_asdomain.notes)
                        else:
                            self.assertEqual(value, getattr(mod_asdomain, key))


    def test_PFAM_domain(self):
        """Check if all the qualifiers are properly stored in PFAM_domain"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_pfams = [i for i in bp_rec.features if i.type == 'PFAM_domain']
        mod_pfams = rec.get_PFAM_domains()
        qualifiers_as_list = ['note', 'label', 'db_xref']
        for bp_pfam, mod_pfam in zip(bp_pfams, mod_pfams):
            for key, value in bp_pfam.qualifiers.items():
                if value is not None and value:
                    #label and db_xref are lists
                    if key not in qualifiers_as_list:
                        if not hasattr(mod_pfam, key):
                            raise AttributeError('%s is not a member of PFAM_domain'%key)
                        self.assertEqual(str(value[0]), str(getattr(mod_pfam, key)))
                    else:
                        if key == 'note':
                            #note is modified to notes in secmet
                            self.assertEqual(value, mod_pfam.notes)
                        else:
                            self.assertEqual(value, getattr(mod_pfam, key))
