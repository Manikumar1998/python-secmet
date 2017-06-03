# vim :set et sts=4 sw=4 fileencoding=utf-8 :
# Licensed under the APL2, see LICENSE for details
"""Secondary Metabolite Record Objects"""

from Bio import SeqIO


class Record(object):
      """A record containing secondary metabolite clusters"""

      def __init__(self, seq_record=None):
            """Initialise a secondary metabolite record

            :param seq_record:  :class:`Bio.SeqRecord.SeqRecord` to read
            :type seq_record:   :class:`Bio.SeqRecord.SeqRecord`
            """
            self._record = seq_record


      @classmethod
      def from_genbank(cls, filename):
            """Initialise a record from a GenBank file

            :param string filename:    file name of the GenBank file to read
            """
            seq_record = SeqIO.read(filename, 'genbank')
            rec = cls(seq_record=seq_record)
            return rec

      @classmethod
      def from_file(cls, filename,filetype):

            """Initialise a record from a file of specified type

            :param string filename:    file name of the file to read
            :param string filetype:    Type of the inputfile
            """
            
            filetype_list = ['gb','genbank','fasta','fas','fa','emb','embl']
            if filetype in filetype_list:
                        if filetype == 'gb' or filetype == 'genbank':
                                type_of_file = 'genbank'
                        elif filetype == 'fas' or filetype == 'fa' or filetype =='fasta':
                                type_of_file = 'fasta'
                        else:
                                type_of_file = 'embl'
                                
                        seq_record = SeqIO.read(filename, type_of_file)
                        rec = cls(seq_record=seq_record)
                        return rec
            else:
                      return None

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
      def clusters(self):
            """A list of secondary metabolite clusters present in the record"""
            if self._record is None:
                return []
            clusters = [i for i in self._record.features if i.type == 'cluster']
            return clusters

      @property
      def gene(self):
            """A list of secondary metabolite genes present in the record"""
            if self._record is None:
                  return []
            gene_list =[i for i in self._record.features if i.type == 'gene']
            return gene_list

      @property
      def CDS(self):
            """A list of secondary metabolite CDS present in the record"""
            if self._record is None:
                  return []
            CDS = [i for i in self._record.features if i.type == 'CDS']
            return CDS
      @property
      def CDS_motif(self):
                """A list of secondary metabolite cds_motifs present in the record"""
                if self._record is None:
                                return []
                cds_motifs_list =[i for i in self._record.features if i.type == 'CDS_motif']
                return cds_motifs_list

      @property
      def source(self):
                if self._record is not None:
                        for i in self._record.features:
                                if i.type == 'source':
                                       return i
                else:
                        return None

      @property
      def feature_types(self):
                type_features =[]
                for i in self._record.features:
                        if i.type not in type_features:
                                type_features.append(i.type)
                return type_features

      def get_cds_from_gene(self,gene_list):
            """Give the CDS corresponding to a particular gene"""
            cds_list =[]
            for gene in gene_list:
                  if type(gene) != type(self.gene[0]):
                        return None
                  else:
                        gene_name = gene.qualifiers.__getattribute__.__self__['gene'][0]
                        cds = self.CDS
                        for i in cds:
                              if i.qualifiers.__getattribute__.__self__['gene'][0] == gene_name:
                                    cds_list.append(i)
            return cds_list


