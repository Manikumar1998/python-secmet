# vim :set et sts=4 sw=4 fileencoding=utf-8 :
# Licensed under the APL2, see LICENSE for details
"""Secondary Metabolite Record Objects"""

import Bio
from Bio import SeqIO


class Record(object):
    """A record containing secondary metabolite clusters"""

    def __init__(self, seq_record=None):
        """Initialise a secondary metabolite record

        :param seq_record:  :class:`Bio.SeqRecord.SeqRecord` to read
        :type seq_record:   :class:`Bio.SeqRecord.SeqRecord`
        """
        self._record = seq_record
        self._cluster_cds = {}  #Dictionary to create cluster-cds hierarchy

    @classmethod
    def from_file(cls, filename, filetype):

        """Initialise a record from a file of specified type

        :param string filename:    file name of the file to read
        :param string filetype:    Type of the inputfile
        """
        if filetype in ['gb', 'genbank']:
            type_of_file = 'genbank'
        elif filetype in ['fa', 'fas', 'fasta']:
            type_of_file = 'fasta'
        elif filetype in ['emb', 'embl']:
            type_of_file = 'embl'
        else:
            return None
        seq_record = SeqIO.read(filename, type_of_file)
        rec = cls(seq_record=seq_record)
        return rec

    @property
    def id(self):
        """Pass through to seq_record object if available"""
        if self._record is not None:
            return self._record.id
        else:
            return "NO_ID_ASSIGNED"

    @property
    def seq(self):
        """Pass through to seq_record object if available"""
        if self._record is not None:
            return self._record.seq
        else:
            return None

    @property
    def annotations(self):
        """Pass through to seq_record object if available"""
        if self._record is not None:
            return self._record.annotations
        else:
            return {}

    @property
    def description(self):
        """Pass through to seq_record object if available"""
        if self._record is not None:
            return self._record.description
        else:
            return ""

    @property
    def feature_types(self):
        """Returns a list of all types of features present in the record"""
        type_features = []
        for i in self._record.features:
            if i.type not in type_features:
                type_features.append(i.type)
        return type_features

    def _features_by_type(self, f_type):
        """Returns a list of features of specified f_type in the record
        param string f_type: Name of the feature
        """
        if f_type in self.feature_types:
            return [i for i in self._record.features if i.type == f_type]
        else:
            return []

    def get_cds_from_gene(self, genes):
        """Returns a list of CDS features corresponding to a list of gene features
            :param list genes : List of gene features
        """
        cds_list = []
        cdss = self._features_by_type('CDS')
        for gene in genes:
            if not isinstance(gene, Bio.SeqFeature.SeqFeature):
                return None
            else:
                gene_name = gene.qualifiers['gene'][0]
                for cds in cdss:
                    if cds.qualifiers['gene'][0] == gene_name:
                        cds_list.append(cds)
                        cdss.remove(cds)    #Removing to reduce the number of operations
        return cds_list

    def get_gene_from_cds(self, cdss):
        """Returns a list of gene features corresponding to a list of CDS features
            :param cdss : List of CDS features
        """
        gene_list = []
        genes = self._features_by_type('gene')
        for cds in cdss:
            if not isinstance(cds, Bio.SeqFeature.SeqFeature):
                return None
            else:
                cds_name = cds.qualifiers['product'][0]
                for gene in genes:
                    if gene.qualifiers['gene'][0].lower() == cds_name.lower():
                        gene_list.append(gene)
                        genes.remove(gene)    #Removing to reduce the number of operations
        return gene_list

    def make_cluster_cds_pair(self, cluster_object, cds_list):
        """Creates a dictionary of cluster objects with corresponding cds objects
              :param cluster_object:  A cluster feature object
              :param cds_list :  list of cds objects corresponding to the cluster object
        """
        self._cluster_cds[cluster_object] = cds_list

    def get_cds_from_cluster(self, cluster_object):
        """Returns the list of CDS feature objects of the given cluster_object
            :param cluster_object: cluster feature object
        """
        return self._cluster_cds[cluster_object]

    def get_cluster_from_cds(self, cds_object):
        """Returns the cluster feature object of the given cds_object
          :param CDS object: CDS feature object
        """
        for cluster, cds_list in self._cluster_cds.items():
            if cds_object in cds_list:
                return cluster
