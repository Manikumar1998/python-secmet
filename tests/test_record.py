from os import path
import Bio
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from secmet.record import Record
from secmet.record import GenericFeature, ClusterFeature, CDSFeature

#Global variables for test file name and its type
filename = 'nisin.gbk'
filetype = 'genbank'

def get_testfile():
    """File path for testing"""
    return path.join(path.dirname(__file__), 'data', filename)

def test_from_file():
    """Test file operations in Record"""
    testfile = get_testfile()
    bp_rec = SeqIO.read(testfile, filetype)
    rec = Record.from_file(testfile)
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
    rec = Record.from_file(testfile)
    assert isinstance(rec.from_biopython(rec._record), Record)

def test_to_biopython():
    """Test to_biopython() in Record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile)
    assert isinstance(rec.to_biopython(), Bio.SeqRecord.SeqRecord)

def test_get_clusters():
    """Test get_clusters() in Record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile)
    bp_rec = SeqIO.read(testfile, filetype)
    bp_clusters = [i for i in bp_rec.features if i.type == 'cluster']
    mod_clusters = [i.to_biopython()[0] for i in rec.get_clusters()]
    assert len(mod_clusters) == len(bp_clusters)
    for bcluster, mcluster in zip(bp_clusters, mod_clusters):
        assert isinstance(mcluster, Bio.SeqFeature.SeqFeature)
        assert bcluster.type == mcluster.type
        assert bcluster.location.__str__() == mcluster.location.__str__()
        for key, value in bcluster.qualifiers.items():
            assert value == mcluster.qualifiers[key]

def test_get_CDSs():
    """Test get_CDSs() in Record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile)
    bp_rec = SeqIO.read(testfile, filetype)
    bp_CDSs = [i for i in bp_rec.features if i.type == 'CDS']
    mod_CDSs = [i.to_biopython()[0] for i in rec.get_CDSs()]
    assert len(mod_CDSs) == len(bp_CDSs)
    for bcds, mcds in zip(bp_CDSs, mod_CDSs):
        assert isinstance(mcds, Bio.SeqFeature.SeqFeature)
        assert bcds.type == mcds.type
        assert bcds.location.__str__() == mcds.location.__str__()
        for key, value in bcds.qualifiers.items():
            assert value == mcds.qualifiers[key]

def test_modified_generic():
    """Test _modified_generic list in Record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile)
    bp_rec = SeqIO.read(testfile, filetype)
    bp_gens = [i for i in bp_rec.features if i.type != 'CDS' and i.type != 'cluster']
    mod_gens = [i.to_biopython()[0] for i in rec._modified_generic]
    assert len(mod_gens) == len(bp_gens)
    for bgen, mgen in zip(bp_gens, mod_gens):
        assert isinstance(mgen, Bio.SeqFeature.SeqFeature)
        assert bgen.type == mgen.type
        assert bgen.location.__str__() == mgen.location.__str__()
        for key, value in bgen.qualifiers.items():
            assert value == mgen.qualifiers[key]

def test_get_cluster_number():
    """Test get_cluster_number() in Record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile)
    clusters = rec.get_clusters()
    for index, cluster in enumerate(clusters):
        assert rec.get_cluster_number(cluster) == index+1

def test_add_feature():
    """Test add_feature() in Record"""
    testfile = get_testfile()
    rec = Record.from_file(testfile)
    no_of_clusters = len(rec.get_clusters())
    no_of_cdss = len(rec.get_CDSs())
    no_of_generics = len(rec._modified_generic)
    new_cluster = ClusterFeature()

    #ClusterFeature should have valid location for adding
    new_cluster.location = FeatureLocation(15100, 15200)
    new_cds = CDSFeature()
    new_generic = GenericFeature()
    rec.add_feature(new_cluster)
    rec.add_feature(new_cds)
    rec.add_feature(new_generic)
    clusters = rec.get_clusters()
    assert no_of_clusters+1 == len(clusters)
    assert no_of_cdss+1 == len(rec.get_CDSs())
    assert no_of_generics+1 == len(rec._modified_generic)
    for index, cluster in enumerate(clusters):
        assert cluster.get_cluster_number() == index+1
