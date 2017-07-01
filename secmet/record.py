# vim :set et sts=4 sw=4 fileencoding=utf-8 :
# Licensed under the APL2, see LICENSE for details
"""Secondary Metabolite Record Objects"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

def cmp_feature_location(a, b):
    "Compare two features by their start/end locations"
    ret = cmp(a.location.start, b.location.start)
    if ret != 0:
        return ret
    return cmp(a.location.end, b.location.end)

def sort_features(seq_record):
    "Sort features in a seq_record by their position"
    #Check if all features have a proper location assigned
    for feature in seq_record.features:
        assert feature.location is not None
    #Sort features by location
    seq_record.features.sort(cmp=cmp_feature_location)

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
        if not isinstance(category, str) and isinstance(info, str):
            raise TypeError("Type of qualifiers should be 'str'")
        if category not in self._qualifiers:
            self._qualifiers[category] = [info]
        else:
            self._qualifiers[category].append(info)
        return None

    def get_qualifier(self, category):
        """Returns a qualifier of given category"""
        if category in self._qualifiers:
            return self._qualifiers[category]

    def to_biopython(self):
        """Returns a Bio.SeqFeature.SeqFeature of given type of feature"""
        if not isinstance(self.location, FeatureLocation):
            raise ValueError("location should be an instance of Bio.SeqFeature.FeatureLocation")
        new_Generic = SeqFeature(FeatureLocation(self.location.start, self.location.end, self.location.strand), type=self.type)
        new_Generic.qualifiers = self._qualifiers.copy()
        return [new_Generic]


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
        self.translation = None
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

            if 'translation' in self._qualifiers:
                self.translation = self._qualifiers['translation'][0]
            self.location = feature.location

    def get_id(self):
        """Returns the id of the CDSFeature"""
        return self.gene

    def get_cluster(self):
        """Returns a ClusterFeature"""
        return self.cluster

    def to_biopython(self):
        """Returns a Bio.SeqFeature.SeqFeature object with all its members"""
        if not isinstance(self.location, FeatureLocation):
            raise ValueError("location should be an instance of Bio.SeqFeature.FeatureLocation")
        self._qualifiers['sec_met'] = self.sec_met
        self._qualifiers['locus_tag'] = [str(self.locus_tag)]
        self._qualifiers['product'] = [str(self.product)]
        self._qualifiers['protein_id'] = [str(self.protein_id)]
        self._qualifiers['gene'] = [str(self.gene)]
        self._qualifiers['translation'] = [str(self.translation)]
        new_CDS = SeqFeature(FeatureLocation(self.location.start, self.location.end, self.location.strand), type=self.type)
        new_CDS.qualifiers = self._qualifiers.copy()
        return [new_CDS]


class ClusterFeature(Feature):
    """A ClusterFeature which subclasses Feature"""
    def __init__(self, feature=None):
        """Initialise a ClusterFeature
            param feature: class 'Bio.SeqFeature.SeqFeature'
        """
        super(ClusterFeature, self).__init__()
        self.contig_edge = None
        self.detection = None
        self.products = []
        self._qualifiers = {}
        self.parent_record = None
        self.type = 'cluster'
        self.note = []

        if feature is not None:
            self._qualifiers = feature.qualifiers
            if 'cutoff' in self._qualifiers:
                self.cutoff = int(self._qualifiers['cutoff'][0])

            if 'extension' in self._qualifiers:
                self.extension = int(self._qualifiers['extension'][0])

            if 'contig_edge' in self._qualifiers:
                self.contig_edge = self._qualifiers['contig_edge'][0]

            if 'note' in self._qualifiers:
                note_list = self._qualifiers['note']
                self.clusternumber = int(note_list[0].split(':')[1])
                if len(note_list) > 1:
                    self.detection = note_list[1]
                if len(note_list) > 2:
                    for i in range(2, len(note_list)):
                        self.note.append(note_list[i])

            if 'product' in self._qualifiers:
                self.products = self._qualifiers['product']
            self.location = feature.location

        self.cdss = []  #At present they are manually assigned for checking

    def _get_cutoff(self):
        try:
            return self.__cutoff
        except:
            return None
    def _set_cutoff(self, value):
        if not isinstance(value, int):
            raise TypeError("cutoff must be an integer")
        self.__cutoff = value
    cutoff = property(_get_cutoff, _set_cutoff)

    def _get_extension(self):
        try:
            return self.__extension
        except:
            return None
    def _set_extension(self, value):
        if not isinstance(value, int):
            raise TypeError("extension must be an integer")
        self.__extension = value
    extension = property(_get_extension, _set_extension)

    def add_product(self, product_value):
        """Adds a product qualifier to the ClusterFeature object"""
        if not isinstance(product_value, str):
            raise TypeError("Type of products should be 'str'")
        self.products.append(product_value)

    def get_products(self):
        """Returns product qualifier from ClusterFeature object"""
        return self.products

    def get_cluster_number(self):
        """Returns the clusternumber of the cluster"""
        if self.parent_record is None:
            raise ValueError('Parent record is None')
        return self.parent_record.get_cluster_number(self)

    def get_CDSs(self):
        """Retruns a list of CDS objects which belong to this cluster"""
        return self.cdss

    def to_biopython(self):
        """Returns a Bio.SeqFeature.SeqFeature object with all its members"""
        if not isinstance(self.location, FeatureLocation):
            raise ValueError("location should be an instance of Bio.SeqFeature.FeatureLocation")
        self._qualifiers['note'] = ["Cluster number: " + str(self.get_cluster_number())]
        self._qualifiers['note'].append(self.detection)
        self._qualifiers['note'].extend(self.note)
        self._qualifiers['cutoff'] = [str(self.cutoff)]
        self._qualifiers['extension'] = [str(self.extension)]
        self._qualifiers['product'] = self.products
        self._qualifiers['contig_edge'] = [str(self.contig_edge)]
        new_Cluster = SeqFeature(FeatureLocation(self.location.start, self.location.end, self.location.strand), type=self.type)
        new_Cluster.qualifiers = self._qualifiers.copy()
        return [new_Cluster]


class Record(object):
    """A record containing secondary metabolite clusters"""

    def __init__(self, seq_record=None):
        """Initialise a secondary metabolite record

        :param seq_record:  :class:`Bio.SeqRecord.SeqRecord` to read
        :type seq_record:   :class:`Bio.SeqRecord.SeqRecord`
        """
        self._record = seq_record
        self._modified_cds = []        #A list containing instances of CDSFeature
        self._modified_cluster = []    #A list containing instances of ClusterFeature
        self._modified_generic = []    #A list containing instances of GenericFeature
        self._cluster_number_dict = {} #A dictionary to map clusters and their numbers

        if self._record is not None:
            if not isinstance(self._record, SeqRecord):
                raise ValueError("SeqRecord should be an instance of 'Bio.SeqRecord.SeqRecord'")
            self.from_biopython(self._record)
        else:
            self._record = SeqRecord(Seq(""))

    @classmethod
    def from_file(cls, filename):

        """Initialise a record from a file of specified type

        :param string filename:    file name of the file to read
        :param string filetype:    Type of the inputfile
        """
        filetype = filename.split('.')[-1]
        if filetype in ['gb', 'gbk', 'genbank']:
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
    @id.setter
    def id(self, value):
        """Setter for id in Record"""
        if not isinstance(value, str):
            raise ValueError('ID should be of type "str"')
        self._record.id = value

    @property
    def seq(self):
        """Pass through to seq_record object if available"""
        if self._record is not None:
            return self._record.seq
        else:
            return None
    @seq.setter
    def seq(self, value):
        """Setter for seq in Record"""
        if not isinstance(value, Seq):
            raise ValueError('Sequence should be of type "Bio.Seq.Seq"')
        self._record.seq = value

    @property
    def description(self):
        """Pass through to seq_record object if available"""
        if self._record is not None:
            return self._record.description
        else:
            return ""
    @description.setter
    def description(self, value):
        """Setter for description in Record"""
        if not isinstance(value, str):
            raise ValueError('Description should be of type "string"')
        self._record.description = value

    @property
    def name(self):
        """Pass through to seq_record object if available"""
        if self._record is not None:
            return self._record.name
        else:
            return "NO_NAME_ASSIGNED"
    @name.setter
    def name(self, value):
        """Setter for name in Record"""
        if not isinstance(value, str):
            raise ValueError('Name should be of type "string"')
        self._record.name = value

    @property
    def annotations(self):
        """Pass through to seq_record object if available"""
        if self._record is not None:
            return self._record.annotations
        else:
            return {}
    def add_annotation(self, key, value):
        """Adding annotations in Record"""
        if not (isinstance(key, str) and (isinstance(value, str) or isinstance(value, list))):
            raise ValueError('Key and Value are not in right format')
        self._record.annotations[key] = value


    def get_clusters(self):
        """A list of secondary metabolite clusters present in the record"""
        return self._modified_cluster
    def set_clusters(self, clusters_list):
        """To set the clusters of the seq_record"""
        self._modified_cluster = clusters_list

    def get_CDSs(self):
        """A list of secondary metabolite clusters present in the record"""
        return self._modified_cds

    def get_generics(self):
        """A list of secondary metabolite generics present in the record"""
        return self._modified_generic
    def set_generics(self, generics_list):
        """To set the generic features of the seq_record"""
        self._modified_generic = generics_list

    def to_biopython(self):
        """Returns a Bio.SeqRecord instance of the record"""
        new_record = self._record
        features = self.get_generics()[:] #Clone the private list
        features.extend(self.get_clusters())
        features.extend(self.get_CDSs())
        record_features = []
        for feature in features:
            record_features.append(feature.to_biopython()[0])
        new_record.features = record_features  #A new_record with all the modified features
        sort_features(new_record)
        return new_record

    def get_cluster_number(self, clusterfeature):
        """Returns cluster number of a cluster feature
            param ClusterFeature clusterfeature : A instance of ClusterFeature class
        """
        return self._cluster_number_dict[clusterfeature]

    def add_feature(self, feature):
        """Adds features to appropriate lists"""
        if not isinstance(feature, Feature):
            raise TypeError("The argument is not an instance of 'Feature'")
        if feature.type == 'cluster':
            if not isinstance(feature.location, FeatureLocation):
                raise ValueError("location should be an instance of Bio.SeqFeature.FeatureLocation")
            clusters = self.get_clusters()
            clusters.append(None)
            for index, cluster in enumerate(clusters):
                if cluster is not None:
                    if feature.location.start < cluster.location.start:
                        break
                else:
                    clusters[index] = feature
                    feature.parent_record = self
                    for index, cluster in enumerate(clusters):
                        self._cluster_number_dict[cluster] = index+1
                    return
            clusters.insert(index, feature)
            feature.parent_record = self
            for index, cluster in enumerate(clusters):
                self._cluster_number_dict[cluster] = index+1
            return

        elif feature.type == 'CDS':
            self._modified_cds.append(feature)
        else:
            self._modified_generic.append(feature)

    def from_biopython(self, record):
        """Modifies _modified_features list with new Feature instances"""
        features = record.features
        for feature in features:
            if feature.type == 'CDS':
                feature = CDSFeature(feature)
                self._modified_cds.append(feature)
            elif feature.type == 'cluster':
                feature = ClusterFeature(feature)
                feature.parent_record = self
                self._modified_cluster.append(feature)
                self._cluster_number_dict[feature] = self._modified_cluster.index(feature)+1
            else:
                feature = GenericFeature(feature)
                self._modified_generic.append(feature)
        return self
