from os import path
import Bio
from Bio import SeqIO
from secmet.record import Record
from secmet.record import GenericFeature, ClusterFeature, CDSFeature


def get_testfile(filename):
    return path.join(path.dirname(__file__), 'data', filename)


def test_init_empty():
    rec = Record()
    assert isinstance(rec, Record)
    assert rec.id == "NO_ID_ASSIGNED"
    assert rec.seq is None
    assert rec.annotations == {}
    assert rec.description == ""
    assert rec.clusters == []

def test_from_file():
    testfile = get_testfile('nisin.gbk')
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
	testfile = get_testfile('nisin.gbk')
	rec = Record.from_file(testfile, 'genbank')
	assert isinstance(rec.from_biopython(), Record)

def test_to_biopython():
	testfile = get_testfile('nisin.gbk')
	rec = Record.from_file(testfile, 'genbank')
	assert isinstance(rec.to_biopython, Bio.SeqRecord.SeqRecord)

def test_get_clusters():
    testfile = get_testfile('nisin.gbk')
    rec = Record.from_file(testfile, 'genbank')
    bp_rec = SeqIO.read(testfile, 'genbank')
    clusters = [i for i in bp_rec.features if i.type == 'cluster']
	#Should call from_biopython() to access features
    rec.from_biopython()
    assert len(rec.get_clusters) == len(clusters)
    assert isinstance(rec.get_clusters[0], ClusterFeature)
    
def test_get_CDSs():
    testfile = get_testfile('nisin.gbk')
    rec = Record.from_file(testfile, 'genbank')
    bp_rec = SeqIO.read(testfile, 'genbank')
    CDSs = [i for i in bp_rec.features if i.type == 'CDS']
	#Should call from_biopython() to access features
    rec.from_biopython()
    assert len(rec.get_CDSs) == len(CDSs)
    assert isinstance(rec.get_CDSs[0], CDSFeature)
