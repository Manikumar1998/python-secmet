from os import path
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from secmet.record import Record
from secmet.record import ClusterFeature

filename = 'nisin.gbk'
filetype = 'genbank'

def get_testfile():
    """File path for testing"""
    return path.join(path.dirname(__file__), 'data', filename)

def test_add_new_cluster():
    """Test for adding a new cluster to record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile)
    new_cluster = ClusterFeature()
    new_cluster.location = FeatureLocation(100, 500)
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

def test_add_existing_cluster():
    """Test for accessing the existing cluster from record"""
    testfile = get_testfile()
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

def write_to_genbank_file():
    """Write data from test_add_new_cluster()"""
    testfile = get_testfile()
    rec = Record.from_file(testfile)
    new_cluster_feature = test_add_new_cluster()
    rec.add_feature(new_cluster_feature)
    record_1 = rec.to_biopython()
    with open('test_'+filename, 'w') as handle:
        SeqIO.write([record_1], handle, filetype)

    #Write data from test_add_existing_cluster()
    rec = Record.from_file(testfile)
    try:
        new_cluster_feature = test_add_existing_cluster()
        rec.add_feature(new_cluster_feature)
    except TypeError:   #To return if no clusters are already present in the file
        return
    record_2 = rec.to_biopython()
    with open('test_'+filename, 'w') as handle:
        SeqIO.write([record_2], handle, filetype)
