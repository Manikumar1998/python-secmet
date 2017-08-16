from os import path
import unittest
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from secmet.record import Record, ClusterFeature

filename = 'Y16952.3.final.gbk'
filetype = 'genbank'

class TestClusterFeature(unittest.TestCase):
    def get_testfile(self):
        """File path for testing"""
        return path.join(path.dirname(__file__), 'data', filename)

    def BioFeature(self):
        biofeature = SeqFeature(location=FeatureLocation(10, 100))
        biofeature.qualifiers = {'contig_edge': ['fake_contig_edge'], 'detection': ['fake_detection'],\
                                 'product': ['fake_products'], 'structure': ['fake_structure'],\
                                 'note': ['Cluster number: 1', 'Detection rule(s): fake_detection', 'fake_notes'],\
                                 'probability': ['fake_probability'], 'subclusterblast': ['fake_subclusterblast'],\
                                 'knownclusterblast': ['fake_knownclusterblast'], 'clusterblast': ['fake_clusterblast'],\
                                 'unknown_qualifier': ['fake_qualifier']}
        return biofeature

    def test_ClusterFeature_members(self):
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        bp_rec = SeqIO.read(testfile, filetype)
        bp_clusters = [i for i in bp_rec.features if i.type == 'cluster']
        mod_clusters = rec.get_clusters()
        #Segregate out qualifiers that are stored in list form
        qualifiers_as_list = ['note', 'product', 'clusterblast', 'subclusterblast', \
                              'knownclusterblast']
        for bp_cluster, mod_cluster in zip(bp_clusters, mod_clusters):
            for key, value in bp_cluster.qualifiers.items():
                if value is not None and value:
                    #clusterblast, subclusterblast and knownclusterblast are lists
                    if key not in qualifiers_as_list:
                        if not hasattr(mod_cluster, key):
                            raise AttributeError('%s is not a member of ClusterFeature'%key)
                        self.assertEqual(str(value[0]), str(getattr(mod_cluster, key)))
                    else:
                        if key == 'note':
                            #notes will not contain 'Cluster number: ' and 'Detection rules: ''
                            self.assertEqual(len(value)-2, len(mod_cluster.notes))
                        elif key == 'product':
                            #product is modified to products in secmet
                            self.assertEqual(bp_cluster.qualifiers['product'], mod_cluster.get_products())
                        else:
                            self.assertEqual(value, getattr(mod_cluster, key))
        cluster = ClusterFeature(FeatureLocation(100, 1000))
        #cutoff, extension should be numbers
        try:
            cluster.cutoff = '-a5000'
        except TypeError:
            pass
        try:
            cluster.extension = 'a5000'
        except TypeError:
            pass

        #If valid qualifiers and values are added, We shouldn't get an error
        try:
            cluster.cutoff = 50000
            cluster.extension = 50000
        except:
            raise RuntimeError('Secmet unable to add valid qualifiers')

    def test_BioFeature_to_ClsuterFeature(self):
        biofeature = self.BioFeature()
        cluster_feature = ClusterFeature(feature=biofeature)
        self.assertEqual(str(cluster_feature.location), str(FeatureLocation(10, 100)))
        self.assertEqual(cluster_feature.type, 'cluster')
        self.assertEqual(cluster_feature.contig_edge, 'fake_contig_edge')
        self.assertEqual(cluster_feature.detection, 'Detection rule(s): fake_detection')
        self.assertEqual(cluster_feature.products, ['fake_products'])
        self.assertEqual(cluster_feature.structure, 'fake_structure')
        self.assertEqual(cluster_feature.probability, 'fake_probability')
        self.assertEqual(cluster_feature.subclusterblast, ['fake_subclusterblast'])
        self.assertEqual(cluster_feature.knownclusterblast, ['fake_knownclusterblast'])
        self.assertEqual(cluster_feature.clusterblast, ['fake_clusterblast'])
        self.assertEqual(cluster_feature.notes, ['fake_notes'])
        self.assertEqual(cluster_feature._qualifiers['unknown_qualifier'], ['fake_qualifier'])
        self.assertEqual(repr(cluster_feature), repr(cluster_feature.to_biopython()[0]))
        self.assertIsInstance(cluster_feature.to_biopython()[0], SeqFeature)

    def test_add_product(self):
        cluster = ClusterFeature(FeatureLocation(1000, 10000))
        self.assertEqual([], cluster.get_products())
        try:
            cluster.add_product(111)
        except TypeError:
            pass
        cluster.add_product('fake_product')
        self.assertEqual(['fake_product'], cluster.get_products())

    def test_add_new_cluster(self):
        """Test for adding a new cluster to record"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        new_cluster = ClusterFeature(FeatureLocation(100, 500))
        try:
            new_cluster.cutoff = 300
        except:
            raise ValueError('Error assigning cutoff value')
        try:
            new_cluster.extension = 300
        except:
            raise ValueError('Error assiging extension value')
        new_cluster.contig_edge = True
        new_cluster.detection = 'Detection rules...'
        new_cluster.add_product('product_info')
        no_clusters_initial = len(rec.get_clusters())
        rec.add_feature(new_cluster)
        no_clusters_final = len(rec.get_clusters())
        assert no_clusters_initial+1 == no_clusters_final
        return new_cluster

    def test_add_existing_cluster(self):
        """Test for accessing the existing cluster from record"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        if len(rec.get_clusters()) >= 1:
            new_cluster = rec.get_clusters()[0]
            assert isinstance(new_cluster, ClusterFeature)
            new_cluster.location = FeatureLocation(100, 500)
            try:
                new_cluster.cutoff = 300
            except:
                raise ValueError('Error assigning cutoff value')
            try:
                new_cluster.extension = 300
            except:
                raise ValueError('Error assiging extension value')
            rec.add_feature(new_cluster)
            return new_cluster

    def test_write_to_file(self):
        """Write data from test_add_new_cluster()"""
        testfile = self.get_testfile()
        rec = Record.from_file(testfile)
        new_cluster_feature = self.test_add_new_cluster()
        rec.add_feature(new_cluster_feature)
        record_1 = rec.to_biopython()
        with open('test_new_'+filename, 'w') as handle:
            SeqIO.write([record_1], handle, filetype)

        #Write data from test_add_existing_cluster()
        rec = Record.from_file(testfile)
        try:
            new_cluster_feature = self.test_add_existing_cluster()
            rec.add_feature(new_cluster_feature)
        except TypeError:   #To return if no clusters are already present in the file
            return
        record_2 = rec.to_biopython()
        with open('test_existing_'+filename, 'w') as handle:
            SeqIO.write([record_2], handle, filetype)
