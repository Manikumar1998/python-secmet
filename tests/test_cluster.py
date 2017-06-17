from os import path
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from secmet.record import Record
from secmet.record import ClusterFeature

filename = 'nisin.gbk'

def get_testfile():
    """File path for testing"""
    return path.join(path.dirname(__file__), 'tests/data', filename)

def test_add_new_cluster():
    """Test for adding a new cluster to record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile, 'genbank')
    new_cluster = ClusterFeature()
    new_cluster.location = FeatureLocation(15100, 15500)
    try:
        new_cluster.cutoff = 15300
    except:
        raise ValueError('Error assigning cutoff value')
    try:
        new_cluster.extension = 15300
    except:
        raise ValueError('Error assiging extension value')
    new_cluster.contig_edge = True
    new_cluster.detection = 'Detection rules...'
    new_cluster.add_product('product_info')
    assert len(rec.get_clusters()) == 1
    rec.add_feature(new_cluster)
    assert len(rec.get_clusters()) == 2
    return new_cluster.to_biopython()

def test_add_existing_cluster():
    """Test for accessing the existing cluster from record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile, 'genbank')
    assert len(rec.get_clusters()) == 1
    new_cluster = rec.get_clusters()[0]
    assert isinstance(new_cluster, ClusterFeature)
    new_cluster.location = FeatureLocation(100, 15106)
    try:
        new_cluster.cutoff = 5000
    except:
        raise ValueError('Error assigning cutoff value')
    try:
        new_cluster.extension = 5000
    except:
        raise ValueError('Error assiging extension value')
    rec.add_feature(new_cluster)
    return new_cluster.to_biopython()

def write_to_genbank_file():
    """Write data from test_add_new_cluster()"""
    testfile = get_testfile()
    rec = Record.from_file(testfile, 'genbank')
    record_1 = rec.to_biopython()
    new_cluster_feature = test_add_new_cluster()[0]
    record_1.features.append(new_cluster_feature)
    with open('test_new_cluster.gbk', 'w') as handle:
        SeqIO.write([record_1], handle, "genbank")

    #Write data from test_add_existing_cluster()
    record_2 = rec.to_biopython()
    new_cluster_feature = test_add_existing_cluster()[0]
    record_2.features.append(new_cluster_feature)
    with open('test_existing_cluster.gbk', 'w') as handle:
        SeqIO.write([record_2], handle, "genbank")
