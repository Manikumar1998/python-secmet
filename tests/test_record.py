from os import path
from Bio import SeqIO
from secmet.record import Record


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


def test_from_genbank():
    testfile = get_testfile('nisin.gbk')
    bp_rec = SeqIO.read(testfile, 'genbank')
    rec = Record.from_genbank(testfile)
    assert isinstance(rec, Record)
    assert rec.id == bp_rec.id
    assert rec.seq == bp_rec.seq
    # SNAG: Can't compare Reference objects in Biopython :(
    # So delete them to make the test work.
    del rec.annotations['references']
    del bp_rec.annotations['references']
    assert rec.annotations == bp_rec.annotations
    assert rec.description == bp_rec.description


def test_clusters():
    testfile = get_testfile('nisin.gbk')
    rec = Record.from_genbank(testfile)
    assert len(rec.clusters) == 1

def test_gene():
    testfile = get_testfile('nisin.gbk')
    rec = Record.from_genbank(testfile)
    bp_rec = SeqIO.read(testfile, 'genbank')
    bp_cds = [i for i in bp_rec.features if i.type == 'gene']
    assert len(bp_cds) == len(rec.gene)

def test_cds():
    testfile = get_testfile('nisin.gbk')
    rec = Record.from_genbank(testfile)
    bp_rec = SeqIO.read(testfile, 'genbank')
    bp_cds = [i for i in bp_rec.features if i.type == 'CDS']
    assert len(bp_cds) == len(rec.CDS)


