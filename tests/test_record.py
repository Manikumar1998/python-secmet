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
    bp_gene = [i for i in bp_rec.features if i.type == 'gene']
    assert len(bp_cds) == len(rec.gene)

def test_cds():
    testfile = get_testfile('nisin.gbk')
    rec = Record.from_genbank(testfile)
    bp_rec = SeqIO.read(testfile, 'genbank')
    bp_cds = [i for i in bp_rec.features if i.type == 'CDS']
    assert len(bp_cds) == len(rec.CDS)

def test_get_cds_from_gene():
		testfile = get_testfile('nisin.gbk')
		rec = Record.from_genbank(testfile)
		bp_rec = SeqIO.read(testfile, 'genbank')
		bp_gene = [i for i in bp_rec.features if i.type == 'gene']
		bp_cds = [i for i in bp_rec.features if i.type == 'CDS']
		#get gene name from bp_gene list
		bp_gene_name = bp_gene[0].qualifiers.__getattribute__.__self__['gene'][0]
		
		#get cds name from bp_cds list
		bp_cds_name = bp_cds[0].qualifiers.__getattribute__.__self__['gene'][0] 
		
		#compare bp_gene_name and bp_cds_name
		assert bp_gene_name == bp_cds_name
		
		#compare bp_cds_name and secmet rec cds name
		assert bp_cds_name == rec.get_cds_from_gene(bp_gene[0]).qualifiers.__getattribute__.__self__['gene'][0]








