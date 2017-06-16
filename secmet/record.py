# vim :set et sts=4 sw=4 fileencoding=utf-8 :
# Licensed under the APL2, see LICENSE for details
"""Secondary Metabolite Record Objects"""

from Bio import SeqIO, SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

class Feature(object):
    """A Feature super class that expands to different subclasses"""
    def __init__(self):
        """ Initialise a feature object"""
        self.location = None
        self.type = None


class GenericFeature(Feature):
    """A GenericFeature Feature subclasses Feature
        (Features other than CDSFeature and ClusterFeature)
    """
    def __init__(self, feature=None):
        """Initialise a GenericFeature
            param feature: class 'Bio.SeqFeature.SeqFeature'
        """
        super(GenericFeature, self).__init__()
        self._qualifiers = {}
        if feature is not None:
            self._qualifiers = feature.qualifiers
            self.location = feature.location
            self.type = feature.type

    def add_qualifier(self, category, info):
        """Adds a qualifier to qualifiers dictionary"""
        if isinstance(category, str) and isinstance(info, str):
            if category not in self._qualifiers:
                self._qualifiers[category] = [info]
            else:
                self._qualifiers[category].append(info)
            return None
        else:
            raise TypeError("Type of qualifiers should be 'str'")

    def get_qualifier(self, category):
        """Returns a qualifier of given category"""
        if category in self._qualifiers:
            return self._qualifiers[category]

    def to_biopython(self):
        """Returns a Bio.SeqFeature.SeqFeature of given type of feature"""
        if isinstance(self.location, FeatureLocation):
            new_Generic = SeqFeature(FeatureLocation(self.location.start, self.location.end), type=self.type)
            new_Generic.qualifiers = self._qualifiers.copy()
            return [new_Generic]
        else:
            raise ValueError("location should be an instance of 'Bio.SeqFeature.FeatureLocation'")

class CDSFeature(Feature):
    """A CDSFeature subclasses Feature"""

    def __init__(self, feature=None):
        """Initialise a CDSFeature
            param feature: class 'Bio.SeqFeature.SeqFeature'
        """
        super(CDSFeature, self).__init__()
        self.sec_met = []
        self.locus_tag = None
        self.product = None
        self.protein_id = None
        self.gene = None
        self.cluster = None  #At present we are manually assigning it for checking
        self._qualifiers = {}
        self.type = 'CDS'

        if feature is not None:

            self._qualifiers = feature.qualifiers

            if 'sec_met' in self._qualifiers:
                self.sec_met = self._qualifiers['sec_met']

            if 'locus_tag' in self._qualifiers:
                self.locus_tag = self._qualifiers['locus_tag'][0]

            if 'product' in self._qualifiers:
                self.product = self._qualifiers['product'][0]

            if 'protein_id' in self._qualifiers:
                self.protein_id = self._qualifiers['protein_id'][0]

            if 'gene' in self._qualifiers:
                self.gene = self._qualifiers['gene'][0]
            self.location = feature.location

    def get_id(self):
        """Returns the id of the CDSFeature"""
        return self.gene

    def get_cluster(self):
        """Returns a ClusterFeature"""
        return self.cluster

    def to_biopython(self):
        """Returns a Bio.SeqFeature.SeqFeature object with all its members"""
        if isinstance(self.location, FeatureLocation):
            self._qualifiers['sec_met'] = self.sec_met
            self._qualifiers['locus_tag'] = [str(self.locus_tag)]
            self._qualifiers['product'] = [str(self.product)]
            self._qualifiers['protein_id'] = [str(self.protein_id)]
            self._qualifiers['gene'] = [str(self.gene)]
            new_CDS = SeqFeature(FeatureLocation(self.location.start, self.location.end), type=self.type)
            new_CDS.qualifiers = self._qualifiers.copy()
            return [new_CDS]
        else:
            raise ValueError("location should be an instance of 'Bio.SeqFeature.FeatureLocation'")

class ClusterFeature(Feature):
    """A ClusterFeature which subclasses Feature"""
    def __init__(self, feature=None):
        """Initialise a ClusterFeature
            param feature: class 'Bio.SeqFeature.SeqFeature'
        """
        super(ClusterFeature, self).__init__()
        self.cutoff = None
        self.extension = None
        self.contig_edge = None
        self.detection = None
        self.products = []
        self._qualifiers = {}
        self.parent_record = None
        self.type = 'cluster'

        if feature is not None:
            self._qualifiers = feature.qualifiers
            if 'cutoff' in self._qualifiers:
                self.cutoff = self._qualifiers['cutoff'][0]

            if 'extension' in self._qualifiers:
                self.extension = self._qualifiers['extension'][0]

            if 'contig_edge' in self._qualifiers:
                self.contig_edge = self._qualifiers['contig_edge'][0]

            if 'note' in self._qualifiers:
                self.detection = self._qualifiers['note'][1]
                self.clusternumber = int(self._qualifiers['note'][0].split(':')[1])

            if 'product' in self._qualifiers:
                self.products = self._qualifiers['product']
            self.location = feature.location

        self.cdss = []  #At present they are manually assigned for checking

    def add_product(self, product_value):
        """Adds a product qualifier to the ClusterFeature object"""
        if isinstance(product_value, str):
            self.products.append(product_value)
        else:
            raise TypeError("Type of products should be 'str'")

    def get_products(self):
        """Returns product qualifier from ClusterFeature object"""
        return self.products

    def get_cluster_number(self):
        """Returns the clusternumber of the cluster"""
        if self.parent_record is not None:
            return self.parent_record.get_cluster_number(self)
        else:
            raise ValueError('Parent record is None')

    def get_CDSs(self):
        """Retruns a list of CDS objects which belong to this cluster"""
        return self.cdss

    def to_biopython(self):
        """Returns a Bio.SeqFeature.SeqFeature object with all its members"""
        if isinstance(self.location, FeatureLocation):
            self._qualifiers['note'] = ["Cluster number: " + str(self.get_cluster_number)]
            self._qualifiers['note'].append(self.detection)
            self._qualifiers['cutoff'] = [str(self.cutoff)]
            self._qualifiers['extension'] = [str(self.extension)]
            self._qualifiers['product'] = self.products
            self._qualifiers['contig_edge'] = [str(self.contig_edge)]
            new_Cluster = SeqFeature(FeatureLocation(self.location.start, self.location.end), type=self.type)
            new_Cluster.qualifiers = self._qualifiers.copy()
            return [new_Cluster]
        else:
            raise ValueError("location should be an instance of 'Bio.SeqFeature.FeatureLocation'")

class Record(object):
    """A record containing secondary metabolite clusters"""

    def __init__(self, seq_record=None):
        """Initialise a secondary metabolite record

        :param seq_record:  :class:`Bio.SeqRecord.SeqRecord` to read
        :type seq_record:   :class:`Bio.SeqRecord.SeqRecord`
        """
        self._record = seq_record
        self._modified_CDS = []      #A list containing instances of CDSFeature
        self._modified_Cluster = [] #A list containing instances of ClusterFeature
        self._modified_Generic = []  #A list containing instances of GenericFeature
        self._cluster_number_dict = {}

        if isinstance(self._record, SeqRecord.SeqRecord):
            self.from_biopython(self._record)
        else:
            raise ValueError("SeqRecord should be an instance of 'Bio.SeqRecord.SeqRecord'")

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

    def get_clusters(self):
        """A list of secondary metabolite clusters present in the record"""
        return self._modified_Cluster

    def get_CDSs(self):
        """A list of secondary metabolite clusters present in the record"""
        return self._modified_CDS

    def to_biopython(self):
        """Returns a Bio.SeqRecord instance of the record"""
        new_record = self._record
        features = self._modified_Generic
        features.extend(self._modified_CDS)
        features.extend(self._modified_Cluster)
        record_features = []
        for feature in features:
            record_features.append(feature.to_biopython[0])
        new_record.features = record_features  #A new_record with all the modified features
        return new_record

    def get_cluster_number(self, clusterfeature):
        """Returns cluster number of a cluster feature
            param ClusterFeature clusterfeature : A instance of ClusterFeature class
        """
        return self._modified_Cluster.index(clusterfeature)+1

    def add_feature(self, feature):
        """Adds features to appropriate lists"""
        if isinstance(feature, Feature):
            if feature.type == 'cluster':
                if isinstance(feature.location, FeatureLocation):
                    clusters = self._modified_Cluster
                    clusters.append(None)
                    for index, cluster in enumerate(clusters):
                        if cluster is not None:
                            if feature.location.start < cluster.location.start and feature.location.end < cluster.location.start:
                                break
                        else:
                            clusters[index] = feature
                            feature.parent_record = self
                            return
                    for k in range(len(clusters)-1, index, -1):
                        temp = clusters[k]
                        clusters[k] = clusters[k-1]
                        clusters[k-1] = temp
                    clusters[k-1] = feature
                    feature.parent_record = self
                else:
                    raise ValueError("location should be an instance of 'Bio.SeqFeatures.FeatureLocation'")

            elif feature.type == 'CDS':
                self._modified_CDS.append(feature)
            else:
                self._modified_Generic.append(feature)
        else:
            raise TypeError("The argument is not an instance of 'Feature'")

    def from_biopython(self, record):
        """Modifies _modified_features list with new Feature instances"""
        features = record.features
        for feature in features:
            if feature.type == 'CDS':
                feature = CDSFeature(feature)
                self._modified_CDS.append(feature)
            elif feature.type == 'cluster':
                feature = ClusterFeature(feature)
                feature.parent_record = self
                self._modified_Cluster.append(feature)
            else:
                feature = GenericFeature(feature)
                self._modified_Generic.append(feature)
        return self
