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
            :param cds_array    :List that has the objects of type cluster_cds
            """
            self._record = seq_record
            self.cds_array =[]      #Initialising the cds_array to None, Maximum 100 cluster features

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

      def get_gene_from_cds(self,cds_list):
            """Returns a list of gene features corresponding to a list of CDS features
                :param cds_list : list of CDS features
            """
            gene_list =[]
            for cds in cds_list:
                  if type(cds) != type(self.CDS[0]):
                        return None
                  else:
                        cds_name = cds.qualifiers['product'][0]
                        gene = self.gene
                        for i in gene:
                              if i.qualifiers.__getattribute__.__self__['gene'][0].lower() == cds_name.lower():  #.lower() is used to overcome strings comparison
                                                                                                                  # without considering the case of the alphabets
                                    gene_list.append(i)
            return gene_list

      def make_cluster_cds_pair(self,cluster_object,cds_list):
            """Links cluster objects with corresponding cds objects
                  :param cluster_cobject:  cluster object
                  :param cds_list :  list of cds feature object linked to the cluster object
            """

            self.cds_array.append(cluster_cds(cds_list,cluster_object.id))

      def get_cds_from_cluster(self,cluster_object):
          """ Given a cluster feature object returns the corresponding CDS features list
              :param cluster_object: cluster feature object
          """
          for i in self.cds_array:
                if i .key == cluster_object.id:
                    return i.cds_list
                else:
                    return None

      def get_cluster_from_cds(self,cds_object):
          """ Given a cds feature object returns the corresponding cluster object
              :param CDS object: CDS feature object
          """
          for cluster_cds_object in self.cds_array:
               for cds_obj in cluster_cds_object.cds_list:
                              if cds_obj.qualifiers['product'][0] == cds_object.qualifiers['product'][0]:
                                    for cluster in self.clusters:
                                          if cluster.id == cluster_cds_object.key:
                                                return cluster



class cluster_cds():
      def __init__(self,cds_list=[],key=None):
            self.cds_list = cds_list
            self.key = key


rec = Record.from_file('../tests/data/nisin.gbk','genbank')
rec.make_cluster_cds_pair(rec.clusters[0],rec.CDS)
print rec.get_gene_from_cds(rec.CDS)
