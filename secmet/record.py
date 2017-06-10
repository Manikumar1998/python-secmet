# vim :set et sts=4 sw=4 fileencoding=utf-8 :
# Licensed under the APL2, see LICENSE for details
"""Secondary Metabolite Record Objects"""

from Bio import SeqIO

#pylint: disable=W0312,C0111

class Feature(object):
	"""A Feature super class that expands to different subclasses"""
	def __init__(self, feature):
		""" Initialise a feature object
			param feature: class 'Bio.SeqFeature.SeqFeature'
		"""
		self.feature = feature
		self._qualifiers = self.feature.qualifiers
		self.location = self.feature.location
		self.type = self.feature.type

	def to_biopython(self):
		return [self.feature]


class GenericFeature(Feature):
	"""A GenericFeature Feature subclasses Feature
		(Features other than CDSFeature and ClusterFeature)
	"""
	def __init__(self, feature=None):
		"""Initialise a GenericFeature
			param feature: class 'Bio.SeqFeature.SeqFeature'
		"""
		if feature != None:
			Feature.__init__(self, feature)

	def to_biopython(self):
		"""Returns Bio.SeqFeature.SeqFeature object of the same feature"""
		return [self.feature]

	def add_qualifier(self, category, info):
		"""Adds a qualifier to qualifiers dictionary"""
		self._qualifiers[category] = [info]
		return None

	def get_qualifier(self, category):
		"""Returns a qualifier of given category"""
		return self._qualifiers[category]


class CDSFeature(Feature):
	"""A CDSFeature subclasses Feature"""

	def __init__(self, feature=None):
		"""Initialise a CDSFeature
			param feature: class 'Bio.SeqFeature.SeqFeature'
		"""
		Feature.__init__(self, feature)
		if 'sec_met' in self._qualifiers:
			self.sec_met = self._qualifiers['sec_met']
		else:
			self.sec_met = None

		if 'locus_tag' in self._qualifiers:
			self.locus_tag = self._qualifiers['locus_tag']
		else:
			self.locus_tag = None

		if 'product' in self._qualifiers:
			self.product = self._qualifiers['product'][0]
		else:
			self.product = None

		if 'protein_id' in self._qualifiers:
			self.protein_id = self._qualifiers['protein_id'][0]
		else:
			self.protein_id = None

		if 'gene' in self._qualifiers:
			self.gene = self._qualifiers['gene'][0]
		else:
			self.gene = None

	def get_id(self):
		"""Returns the id of the CDSFeature"""
		return self.gene

	def to_biopython(self):
		"""Returns a Bio.SeqFeature.SeqFeature object of the same feature"""
		return [self.feature]



class ClusterFeature(Feature):
	"""A ClusterFeature which subclasses Feature"""
	def __init__(self, feature=None):
		"""Initialise a ClusterFeature
			param feature: class 'Bio.SeqFeature.SeqFeature'
		"""
		if feature != None:
			Feature.__init__(self, feature)
		if 'cutoff' in self._qualifiers:
			self.cutoff = self._qualifiers['cutoff']
		if 'extension' in self._qualifiers:
			self.extension = self._qualifiers['extension']
		if 'contig_edge' in self._qualifiers:
			self.contig_edge = self._qualifiers['contig_edge']
		if 'detection' in self._qualifiers:
			self.detection = self._qualifiers['note'][1]

	def to_biopython(self):
		"""Returns a Bio.SeqFeature.SeqFeature object of the same feature"""
		return [self.feature]

	def get_products(self):
		"""Returns product qualifier from ClusterFeature object"""
		return self._qualifiers['product']

	def get_cluster_number(self):
		"""Returns the clusternumber of the cluster"""
		return int(self._qualifiers['note'][0].split()[2])


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
