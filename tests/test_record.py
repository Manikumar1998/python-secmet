from os import path
import unittest
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from secmet.record import Record, GenericFeature, ClusterFeature, CDSFeature \
                                , CDS_motifFeature, aSDomain, PFAM_domain

#Global variables for test file name and its type
filename = 'Y16952.3.final.gbk'
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
        if 'references' in rec.annotations:
            del rec.annotations['references']
        if 'references' in bp_rec.annotations:
            del bp_rec.annotations['references']
        self.assertEqual(len(rec), 66669)
        self.assertEqual(rec.annotations, bp_rec.annotations)
        self.assertEqual(rec.description, bp_rec.description)

    def test_empty_Record(self):
        """Test the identifiers of empty Record"""
        rec = Record()
        #seq should be a instance of Bio.Seq.Seq
        with self.assertRaises(ValueError):
            rec.seq = 'FAKE'
        #description, name and id are strings
        with self.assertRaises(ValueError):
            rec.description = 123
        with self.assertRaises(ValueError):
            rec.name = 123
        with self.assertRaises(ValueError):
            rec.id = 123

        rec.id = "fake_id"
        rec.name = 'fake_name'
        rec.seq = Seq("FAKE")
        rec.description = 'fake_description'
        self.assertEqual(rec.id, 'fake_id')
        self.assertEqual(rec.name, 'fake_name')
        self.assertEqual(rec.seq, Seq("FAKE"))
        self.assertEqual(rec.description, 'fake_description')
        self.assertEqual(rec.annotations, {})
        with self.assertRaises(ValueError):
            rec.add_annotation(12, 34)
        rec.add_annotation('fake_key', 'fake_value')
        self.assertEqual(rec.annotations, {'fake_key': 'fake_value'})

    def test_setters(self):
        """Test setters for features lists"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        self.assertNotEqual(rec.get_CDSs(), [])
        self.assertNotEqual(rec.get_clusters(), [])
        self.assertNotEqual(rec.get_PFAM_domains(), [])
        self.assertNotEqual(rec.get_aSDomains(), [])
        self.assertNotEqual(rec.get_generics(), [])

        rec.erase_CDSs()
        rec.erase_clusters()
        rec.erase_generics()
        rec.erase_CDS_motifs()
        rec.erase_PFAM_domains()
        rec.erase_aSDomains()

        self.assertEqual(rec.get_CDSs(), ())
        self.assertEqual(rec.get_clusters(), ())
        self.assertEqual(rec.get_PFAM_domains(), ())
        self.assertEqual(rec.get_aSDomains(), ())
        self.assertEqual(rec.get_generics(), ())

    def test_from_biopython(self):
        """Test from_biopython() in Record"""
        with self.assertRaises(ValueError):
            rec = Record('fake_record')
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
            self.assertEqual(str(bcluster.location), str(mcluster.location))
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
            self.assertEqual(str(bcds.location), str(mcds.location))
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
            self.assertEqual(str(bgen.location), str(mgen.location))
            for key, value in bgen.qualifiers.items():
                self.assertEqual(value, mgen.qualifiers[key])

    def test_get_CDS_motifs(self):
        """Test get_CDS_motifs() in Record"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_CDS_motifs = [i for i in bp_rec.features if i.type == 'CDS_motif']
        mod_CDS_motifs = [i.to_biopython()[0] for i in rec.get_CDS_motifs()]
        self.assertEqual(len(mod_CDS_motifs), len(bp_CDS_motifs))
        for b_motif, m_motif in zip(mod_CDS_motifs, bp_CDS_motifs):
            self.assertIsInstance(m_motif, Bio.SeqFeature.SeqFeature)
            self.assertEqual(b_motif.type, m_motif.type)
            self.assertEqual(str(b_motif.location), str(m_motif.location))
            for key, value in b_motif.qualifiers.items():
                if value:
                    if key in ['score', 'evalue']:
                        self.assertEqual(float(value[0]), float(m_motif.qualifiers[key][0]))
                    else:
                        self.assertEqual(value, m_motif.qualifiers[key])

    def test_get_PFAM_domains(self):
        """Test get_PFAM_domains() in Record"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_PFAM_domains = [i for i in bp_rec.features if i.type == 'PFAM_domain']
        mod_PFAM_domains = [i.to_biopython()[0] for i in rec.get_PFAM_domains()]
        self.assertEqual(len(mod_PFAM_domains), len(bp_PFAM_domains))
        for b_fam, m_fam in zip(mod_PFAM_domains, bp_PFAM_domains):
            self.assertIsInstance(m_fam, Bio.SeqFeature.SeqFeature)
            self.assertEqual(b_fam.type, m_fam.type)
            self.assertEqual(str(b_fam.location), str(m_fam.location))
            for key, value in b_fam.qualifiers.items():
                if value:
                    if key in ['score', 'evalue']:
                        self.assertEqual(float(value[0]), float(m_fam.qualifiers[key][0]))
                    else:
                        self.assertEqual(value, m_fam.qualifiers[key])

    def test_get_aSDomains(self):
        """Test get_aSDomains() in Record"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_aSDomains = [i for i in bp_rec.features if i.type == 'aSDomain']
        mod_aSDomains = [i.to_biopython()[0] for i in rec.get_aSDomains()]
        self.assertEqual(len(mod_aSDomains), len(bp_aSDomains))
        for b_asdomain, m_asdomain in zip(mod_aSDomains, bp_aSDomains):
            self.assertIsInstance(m_asdomain, Bio.SeqFeature.SeqFeature)
            self.assertEqual(b_asdomain.type, m_asdomain.type)
            self.assertEqual(str(b_asdomain.location), str(m_asdomain.location))
            for key, value in b_asdomain.qualifiers.items():
                if value:
                    if key in ['score', 'evalue']:
                        self.assertEqual(float(value[0]), float(m_asdomain.qualifiers[key][0]))
                    else:
                        self.assertEqual(value, m_asdomain.qualifiers[key])

    def test_get_cluster_number(self):
        """Test get_cluster_number() in Record"""
        rec = Record()
        cluster1 = ClusterFeature(FeatureLocation(500, 1500))
        cluster2 = ClusterFeature(FeatureLocation(5000, 6000))
        cluster3 = ClusterFeature(FeatureLocation(2500, 4000))
        rec.add_feature(cluster2)
        self.assertEqual(1, cluster2.get_cluster_number())
        rec.add_feature(cluster1)
        self.assertEqual(1, cluster1.get_cluster_number())
        self.assertEqual(2, cluster2.get_cluster_number())
        rec.add_feature(cluster3)
        self.assertEqual(1, cluster1.get_cluster_number())
        self.assertEqual(3, cluster2.get_cluster_number())
        self.assertEqual(2, cluster3.get_cluster_number())

    def test_cluster_cds_links(self):
        """Test whether cluster(s) and CDS(s) are properly linked"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_clusters = [i for i in bp_rec.features if i.type == 'cluster']
        bp_cdss = [i for i in bp_rec.features if i.type == 'CDS']
        mod_clusters = rec.get_clusters()
        self.assertEqual(len(bp_clusters), len(mod_clusters))
        for bp_cluster, mod_cluster in zip(bp_clusters, mod_clusters):
            bp_clustercdsfeatures = []
            for cds in bp_cdss:
                if bp_cluster.location.start <= cds.location.start <= bp_cluster.location.end or \
                   bp_cluster.location.start <= cds.location.end <= bp_cluster.location.end:
                    bp_clustercdsfeatures.append(cds)
            cluster_cds_features = rec.get_CDSs()+rec.get_clusters()
            for feature in cluster_cds_features:
                rec._update_cluster_cds_links(feature)
            self.assertEqual(len(bp_clustercdsfeatures), len(mod_cluster.get_CDSs()))
            for bp_cds, mod_cds in zip(bp_clustercdsfeatures, mod_cluster.get_CDSs()):
                self.assertEqual(str(bp_cds.location), str(mod_cds.location))
                self.assertEqual(str(mod_cds.get_cluster().location), str(bp_cluster.location), \
                                 str(mod_cluster.location))

    def test_add_feature_cds(self):
        rec = Record()
        cluster1 = ClusterFeature(FeatureLocation(1, 1000))
        cluster2 = ClusterFeature(FeatureLocation(2000, 3000))
        cluster3 = ClusterFeature(FeatureLocation(4000, 5000))
        cds = CDSFeature(FeatureLocation(4500, 4600))
        rec.add_feature(cds)
        #If no clusters are present None should be returned
        self.assertEqual(cds.get_cluster(), None)
        rec.erase_CDSs()
        rec.add_feature(cluster1)
        rec.add_feature(cluster2)
        rec.add_feature(cluster3)
        rec.add_feature(cds)
        self.assertEqual(cds.get_cluster(), cluster3)

    def test_add_feature(self):
        """Test add_feature() in Record"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        no_of_clusters = len(rec.get_clusters())
        no_of_cdss = len(rec.get_CDSs())
        no_of_generics = len(rec.get_generics())
        no_of_cds_motifs = len(rec.get_CDS_motifs())
        no_of_pfam_domains = len(rec.get_PFAM_domains())
        no_of_asdomains = len(rec.get_aSDomains())
        #Create new Feature's with fake identity and fake location
        invalid_feature = 'INVALID_FEATURE'
        new_cluster = ClusterFeature(FeatureLocation(1000, 2000))
        new_cds = CDSFeature(FeatureLocation(1500, 1700))
        new_generic = GenericFeature(FeatureLocation(350, 450), 'FAKE')
        new_cds_motif = CDS_motifFeature(FeatureLocation(150, 200))
        new_pfam_domain = PFAM_domain(FeatureLocation(500, 600))
        new_asdomain = aSDomain(FeatureLocation(600, 700))
        rec.add_feature(new_cluster)
        rec.add_feature(new_cds)
        rec.add_feature(new_generic)
        rec.add_feature(new_cds_motif)
        rec.add_feature(new_pfam_domain)
        rec.add_feature(new_asdomain)
        clusters = rec.get_clusters()
        with self.assertRaises(TypeError):
            rec.add_feature(invalid_feature)
        self.assertEqual(no_of_clusters+1, len(clusters))
        self.assertEqual(no_of_cdss+1, len(rec.get_CDSs()))
        self.assertEqual(no_of_generics+1, len(rec.get_generics()))
        self.assertEqual(no_of_cds_motifs+1, len(rec.get_CDS_motifs()))
        self.assertEqual(no_of_pfam_domains+1, len(rec.get_PFAM_domains()))
        self.assertEqual(no_of_asdomains+1, len(rec.get_aSDomains()))
        for index, cluster in enumerate(clusters):
            self.assertEqual(cluster.get_cluster_number(), index+1)
