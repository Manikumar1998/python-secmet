# vim :set et sts=4 sw=4 fileencoding=utf-8 :
# Licensed under the APL2, see LICENSE for details
"""Secondary Metabolite Record Objects"""

from Bio import SeqIO

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

	@property
	def to_biopython(self):
		"""Returns a Bio.SeqFeature.SeqFeature object of same feature"""
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

	@property
	def get_id(self):
		"""Returns the id of the CDSFeature"""
		return self.gene

	@property
	def get_cluster(self):
		"""Returns a ClusterFeature"""
		#TO-DO: Should return the corresponding ClusterFeature
		return



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
		else:
			self.cutoff = None
		if 'extension' in self._qualifiers:
			self.extension = self._qualifiers['extension']
		else:
			self.extension = None
		if 'contig_edge' in self._qualifiers:
			self.contig_edge = self._qualifiers['contig_edge']
		else:
			self.contig_edge = None
		if 'detection' in self._qualifiers:
			self.detection = self._qualifiers['note'][1]
		else:
			self.detection = None

	@property
	def get_products(self):
		"""Returns product qualifier from ClusterFeature object"""
		return self._qualifiers['product']

	@property
	def get_cluster_number(self):
		"""Returns the clusternumber of the cluster"""
		return int(self._qualifiers['note'][0].split()[2])

	@property
	def get_CDSs(self):
		#TO-DO: Should return a list of CDSFeatures
		return

	def add_product(self, product):
		"""Adds a product qualifier to the ClusterFeature object"""
		self._qualifiers['product'].append(product)


class Record(object):
	"""A record containing secondary metabolite clusters"""

	def __init__(self, seq_record=None):
		"""Initialise a secondary metabolite record

		:param seq_record:  :class:`Bio.SeqRecord.SeqRecord` to read
		:type seq_record:   :class:`Bio.SeqRecord.SeqRecord`
		"""
		self._record = seq_record
		self._modified_features = [] #A list containing instances of Feature or its subclasses


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
	def get_clusters(self):
		"""A list of secondary metabolite clusters present in the record"""
		if self._record is None:
		    return []

		clusters = [i for i in self._modified_features if i.type == 'cluster']
		return clusters

	@property
	def get_CDSs(self):
		"""A list of secondary metabolite clusters present in the record"""
		if self._record is None:
		    return []

		cdss = [i for i in self._modified_features if i.type == 'CDS']
		return cdss

	@property
	def to_biopython(self):
		"""Returns a Bio.SeqRecord instance"""
		return self._record

	def get_cluster_number(self, clusterfeature):
		"""Returns cluster number of a cluster feature
			param ClusterFeature clusterfeature : A instance of ClusterFeature class
		"""
		return clusterfeature.get_cluster_number

	def from_biopython(self):
		"""Modifies _modified_features list with new Feature instances"""
		features = self._record.features
		for feature in features:
			if feature.type == 'CDS':
				feature = CDSFeature(feature)
				self._modified_features.append(feature)
			elif feature.type == 'cluster':
				feature = ClusterFeature(feature)
				self._modified_features.append(feature)
			else:
				feature = GenericFeature(feature)
				self._modified_features.append(feature)
		return self
