from os import path
from Bio import SeqIO
from secmet.record import Record

def get_testfile(filename):
    return path.join(path.dirname(__file__), 'data', filename)

def test_init_empty():
    rec = Record()
    assert isinstance(rec, Record)
    assert rec.id == "NO_ID_ASSIGNED"
    assert rec.annotations == {}
    assert rec.description == ""
    assert rec.clusters == []


def test_from_genbank():
    testfile = get_testfile('nisin.gbk')
    bp_rec = SeqIO.read(testfile, 'genbank')
    rec = Record.from_genbank(testfile)
    assert isinstance(rec, Record)
    assert rec.id == bp_rec.id
    assert rec.annotations == bp_rec.annotations
    assert rec.description == bp_rec.description


def test_clusters():
    testfile = get_testfile('nisin.gbk')
    rec = Record.from_genbank(testfile)
    assert len(rec.clusters) == 1
