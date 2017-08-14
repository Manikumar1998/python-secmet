from os import path
import unittest
import Bio
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from secmet.record import Record, GenericFeature, ClusterFeature, CDSFeature

#Global variables for test file name and its type
filename = 'nisin.gbk'
filetype = 'genbank'

class TestRecordMethods(unittest.TestCase):

    def get_testfile(self):
        """File path for testing"""
        return path.join(path.dirname(__file__), 'data', filename)

    def test_from_file(self):
        """Test file operations in Record"""
        testfile = self.get_testfile()
        bp_rec = SeqIO.read(testfile, filetype)
        rec = Record.from_file(testfile)
        assert isinstance(rec, Record)
        self.assertEqual(rec.id, bp_rec.id)
        self.assertEqual(rec.seq, bp_rec.seq)
        # SNAG: Can't compare Reference objects in Biopython :(
        # So delete them to make the test work.
        del rec.annotations['references']
        del bp_rec.annotations['references']
        self.assertEqual(rec.annotations, bp_rec.annotations)
        self.assertEqual(rec.description, bp_rec.description)

    def test_from_biopython(self):
        """Test from_biopython() in Record"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        self.assertIsInstance(rec.from_biopython(rec._record), Record)

    def test_to_biopython(self):
        """Test to_biopython() in Record"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        self.assertIsInstance(rec.to_biopython(), Bio.SeqRecord.SeqRecord)

    def test_get_clusters(self):
        """Test get_clusters() in Record"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_clusters = [i for i in bp_rec.features if i.type == 'cluster']
        mod_clusters = [i.to_biopython()[0] for i in rec.get_clusters()]
        self.assertEqual(len(mod_clusters), len(bp_clusters))
        for bcluster, mcluster in zip(bp_clusters, mod_clusters):
            self.assertIsInstance(mcluster, Bio.SeqFeature.SeqFeature)
            self.assertEqual(bcluster.type, mcluster.type)
            self.assertEqual(bcluster.location.__str__(), mcluster.location.__str__())
            for key, value in bcluster.qualifiers.items():
                self.assertEqual(value, mcluster.qualifiers[key])

    def test_get_CDSs(self):
        """Test get_CDSs() in Record"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_CDSs = [i for i in bp_rec.features if i.type == 'CDS']
        mod_CDSs = [i.to_biopython()[0] for i in rec.get_CDSs()]
        self.assertEqual(len(mod_CDSs), len(bp_CDSs))
        for bcds, mcds in zip(bp_CDSs, mod_CDSs):
            self.assertIsInstance(mcds, Bio.SeqFeature.SeqFeature)
            self.assertEqual(bcds.type, mcds.type)
            self.assertEqual(bcds.location.__str__(), mcds.location.__str__())
            for key, value in bcds.qualifiers.items():
                if key != 'sec_met':
                    self.assertEqual(value, mcds.qualifiers[key])

    def test_get_generics(self):
        """Test get_generics() in Record"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        non_generic_features = ['CDS', 'cluster', 'CDS_motif', 'aSDomain', 'PFAM_domain']
        bp_gens = [i for i in bp_rec.features if i.type not in non_generic_features]
        mod_gens = [i.to_biopython()[0] for i in rec.get_generics()]
        self.assertEqual(len(mod_gens), len(bp_gens))
        for bgen, mgen in zip(bp_gens, mod_gens):
            self.assertIsInstance(mgen, Bio.SeqFeature.SeqFeature)
            self.assertEqual(bgen.type, mgen.type)
            self.assertEqual(bgen.location.__str__(), mgen.location.__str__())
            for key, value in bgen.qualifiers.items():
                self.assertEqual(value, mgen.qualifiers[key])

    def test_get_cluster_number(self):
        """Test get_cluster_number() in Record"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        clusters = rec.get_clusters()
        for index, cluster in enumerate(clusters):
            self.assertEqual(rec.get_cluster_number(cluster), index+1)

    def test_add_feature(self):
        """Test add_feature() in Record"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        no_of_clusters = len(rec.get_clusters())
        no_of_cdss = len(rec.get_CDSs())
        no_of_generics = len(rec.get_generics())
        new_cluster = ClusterFeature(FeatureLocation(15100, 15200))
        new_cds = CDSFeature(FeatureLocation(200, 300))
        new_generic = GenericFeature(FeatureLocation(350, 450), 'FAKE')
        rec.add_feature(new_cluster)
        rec.add_feature(new_cds)
        rec.add_feature(new_generic)
        clusters = rec.get_clusters()
        self.assertEqual(no_of_clusters+1, len(clusters))
        self.assertEqual(no_of_cdss+1, len(rec.get_CDSs()))
        self.assertEqual(no_of_generics+1, len(rec._modified_generic))
        for index, cluster in enumerate(clusters):
            self.assertEqual(cluster.get_cluster_number(), index+1)
