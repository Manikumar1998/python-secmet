from os import path
import Bio
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from secmet.record import Record
from secmet.record import GenericFeature, ClusterFeature, CDSFeature


filename = 'nisin.gbk'

def get_testfile():
    """File path for testing"""
    return path.join(path.dirname(__file__), 'tests/data', filename)

def test_from_file():
    """Test file operations in Record"""
    testfile = get_testfile()
    bp_rec = SeqIO.read(testfile, 'genbank')
    rec = Record.from_file(testfile, 'genbank')
    assert isinstance(rec, Record)
    assert rec.id == bp_rec.id
    assert rec.seq == bp_rec.seq
    # SNAG: Can't compare Reference objects in Biopython :(
    # So delete them to make the test work.
    del rec.annotations['references']
    del bp_rec.annotations['references']
    assert rec.annotations == bp_rec.annotations
    assert rec.description == bp_rec.description

def test_from_biopython():
    """Test from_biopython() in Record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile, 'genbank')
    assert isinstance(rec.from_biopython(rec._record), Record)

def test_to_biopython():
    """Test to_biopython() in Record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile, 'genbank')
    assert isinstance(rec.to_biopython(), Bio.SeqRecord.SeqRecord)

def test_get_clusters():
    """Test get_clusters() in Record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile, 'genbank')
    bp_rec = SeqIO.read(testfile, 'genbank')
    bp_clusters = [i for i in bp_rec.features if i.type == 'cluster']
    mod_clusters = rec.get_clusters()
    assert len(mod_clusters) == len(bp_clusters)
    for cluster in mod_clusters:
        assert isinstance(cluster, ClusterFeature)

def test_get_CDSs():
    """Test get_CDSs() in Record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile, 'genbank')
    bp_rec = SeqIO.read(testfile, 'genbank')
    bp_CDSs = [i for i in bp_rec.features if i.type == 'CDS']
    mod_CDSs = rec.get_CDSs()
    assert len(mod_CDSs) == len(bp_CDSs)
    for cds in mod_CDSs:
        assert isinstance(cds, CDSFeature)

def test_get_cluster_number():
    """Test get_cluster_number() in Record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile, 'genbank')
    cluster = rec.get_clusters()[0]
    assert rec.get_cluster_number(cluster) == 1

def test_add_feature():
    """Test add_feature() in Record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile, 'genbank')
    no_of_clusters = len(rec.get_clusters())
    no_of_cdss = len(rec.get_CDSs())
    no_of_generics = len(rec._modified_Generic)
    new_cluster = ClusterFeature()

    #ClusterFeature should have valid location for adding
    new_cluster.location = FeatureLocation(15100, 15200)
    new_cds = CDSFeature()
    new_generic = GenericFeature()
    rec.add_feature(new_cluster)
    rec.add_feature(new_cds)
    rec.add_feature(new_generic)
    assert no_of_clusters+1 == len(rec.get_clusters())
    assert no_of_cdss+1 == len(rec.get_CDSs())
    assert no_of_generics+1 == len(rec._modified_Generic)
    assert rec.get_clusters()[0].get_cluster_number() == 1
    assert new_cluster.get_cluster_number() == 2
