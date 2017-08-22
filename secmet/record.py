# vim :set et sts=4 sw=4 fileencoding=utf-8 :
# Licensed under the APL2, see LICENSE for details
"""Secondary Metabolite Record Objects"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from bisect import bisect_left

def cmp_feature_location(a, b):
    "Compare two features by their start/end locations"
    ret = cmp(a.location.start, b.location.start)
    if ret != 0:
        return ret
    return cmp(a.location.end, b.location.end)

def sort_features(seq_record):
    "Sort features in a seq_record using their locations"
    seq_record.features.sort(cmp=cmp_feature_location)


def find_new_cluster_pos(clusters, target_cluster):
    """Find appropriate position in array to add new cluster using Bisection method
        param clusters: A list of all existing ClusterFeature(s) in the record
        param target_cluster: An instance of ClusterFeature
    """
    if not clusters:
        return 0
    cluster_start_locations = [cluster.location.start for cluster in clusters]
    return bisect_left(cluster_start_locations, target_cluster.location.start)


def find_cluster_of_new_cds(clusters, new_cds):
    """Find the corresponding cluster feature of a cds feature using Bisection method
        param clusters: A list of all existing ClusterFeature(s) in the record
        param new_cds: An instance of CDSFeature
    """
    if not clusters:
        return
    if new_cds.location.end < clusters[0].location.start or \
       new_cds.location.start > clusters[len(clusters)-1].location.end:
        return
    else:
        cluster_starts = [cluster.location.start for cluster in clusters]
        cluster_ends = [cluster.location.end for cluster in clusters]
        if bisect_left(cluster_starts, new_cds.location.start)-1 == bisect_left(cluster_ends, new_cds.location.start):
            index = bisect_left(cluster_ends, new_cds.location.start)
            clusters[index].cdss.append(new_cds)
            new_cds.cluster = clusters[index]

        if bisect_left(cluster_starts, new_cds.location.end)-1 == bisect_left(cluster_ends, new_cds.location.end):
            index = bisect_left(cluster_ends, new_cds.location.end)
            clusters[index].cdss.append(new_cds)
            new_cds.cluster = clusters[index]


class Feature(object):
    """A Feature super class that extends to different subclasses"""
    def __init__(self):
        """ Initialise a Feature object"""
        self.type = None
        self.notes = []

    #Check for a valid feature location
    def _get_location(self):
        return self.__location
    def _set_location(self, value):
        if not isinstance(value, (FeatureLocation, CompoundLocation)):
            raise TypeError("Location must be an instance of 'FeatureLocation' or 'CompoundLocation'")
        self.__location = value
    location = property(_get_location, _set_location)

    def extract(self, parent_seq):
        """Return Feature's seq from its parent's seq"""
        if not self.location:
            raise ValueError("Location is None. Extracting Failed")
        return self.location.extract(parent_seq)


class GenericFeature(Feature):
    """A GenericFeature Feature subclasses Feature
        (Features other than CDS, cluster, CDS_motif, PFAM_domain and aSDomin)
    """
    def __init__(self, f_location=None, f_type=None, feature=None):
        """Initialise a GenericFeature
            param f_location: class 'Bio.SeqFeature.FeatureLocation/CompoundLocation'
            param feature: class 'Bio.SeqFeature.SeqFeature'
        """
        super(GenericFeature, self).__init__()
        self._qualifiers = {}
        self.locus_tag = None
        self.translation = None
        self.gene = None
        self.name = None
        self.seq = None
        self.description = None
        self.sec_met = []

        if feature:
            """Initialise class members(qualifiers) using SeqFeature object"""
            self._qualifiers = feature.qualifiers
            self.type = feature.type
            self.location = feature.location
            self.locus_tag = self._qualifiers.pop('locus_tag', [None])[0]
            self.gene = self._qualifiers.pop('gene', [None])[0]
            self.translation = self._qualifiers.pop('translation', [None])[0]
            self.name = self._qualifiers.pop('name', [None])[0]
            self.seq = self._qualifiers.pop('seq', [None])[0]
            self.description = self._qualifiers.pop('description', [None])[0]
            self.sec_met = self._qualifiers.pop('sec_met', [])
            self.notes = self._qualifiers.pop('note', [])
        else:
            self.location = f_location
            if not isinstance(f_type, str):
                raise ValueError('Type of the feature should be a string')
            self.type = f_type

    def add_qualifier(self, category, info):
        """Adds a qualifier to qualifiers dictionary"""
        if not isinstance(category, str):
            if not isinstance(info, (str, list)):
                if not isinstance(info, (int, float)):
                    raise TypeError("Qualifier category should be str and value should be str or list or number")
                else:
                    info = str(info)
        if category in ['evalue', 'score', 'probability']:
            try:
                info = float(info)
            except:
                raise ValueError('%s should be a number'% category)
        if hasattr(self, category):
            if isinstance(getattr(self, category), list):
                if isinstance(info, list):
                    getattr(self, category).extend(info)
                else:
                    getattr(self, category).append(info)
            else:
                setattr(self, category, info)
        else:
            if category not in self._qualifiers:
                if isinstance(info, list):
                    self._qualifiers[category] = info
                else:
                    self._qualifiers[category] = [info]
            else:
                self._qualifiers[category].append(info)
        return None

    def get_qualifier(self, category):
        """Returns a qualifier of given category"""
        if category in self._qualifiers:
            return self._qualifiers[category]
        elif category.lower() in self._qualifiers:
            return self._qualifiers[category.lower()]
        elif category.upper() in self._qualifiers:
            return self._qualifiers[category.upper()]
        else:
            if hasattr(self, category):
                if getattr(self, category):
                    return [getattr(self, category)]
            elif hasattr(self, category.lower()):
                if getattr(self, category.lower()):
                    return [getattr(self, category.lower())]
        return []

    def to_biopython(self):
        """Returns a Bio.SeqFeature.SeqFeature of given type of feature"""
        new_Generic = SeqFeature(self.location, type=self.type)
        if self.locus_tag:
            new_Generic.qualifiers['locus_tag'] = [str(self.locus_tag)]
        if self.translation:
            new_Generic.qualifiers['translation'] = [str(self.translation)]
        if self.gene:
            new_Generic.qualifiers['gene'] = [str(self.gene)]
        if self.name:
            new_Generic.qualifiers['name'] = [str(self.name)]
        if self.seq:
            new_Generic.qualifiers['seq'] = [str(self.seq)]
        if self.description:
            new_Generic.qualifiers['description'] = [str(self.description)]
        if self.sec_met:
            new_Generic.qualifiers['sec_met'] = self.sec_met
        if self.notes:
            new_Generic.qualifiers['note'] = self.notes
        for key, value in self._qualifiers.items():
            new_Generic.qualifiers[key] = value
        return [new_Generic]

    def __repr__(self):
        """A string representation of biopython generic features"""
        return repr(self.to_biopython()[0])


class CDSFeature(Feature):
    """A CDSFeature subclasses Feature"""

    def __init__(self, f_location=None, feature=None):
        """Initialise a CDSFeature
            param f_location: class 'Bio.SeqFeature.FeatureLocation/CompoundLocation'
            param feature: class 'Bio.SeqFeature.SeqFeature'
        """
        super(CDSFeature, self).__init__()
        self.id = '<unknown id>'
        self.sec_met = SecMetQualifier()
        self.locus_tag = None
        self.product = None
        self.protein_id = None
        self.gene = None
        self.translation = None
        self.cluster = None
        self.transl_table = None
        self.source = None
        self.db_xref = []
        self.EC_number = []
        self.aSProdPred = []
        self.aSASF_choice = []
        self.aSASF_note = []
        self.aSASF_prediction = []
        self.aSASF_scaffold = []

        self._qualifiers = {}
        self.sec_met_predictions = []
        self.type = 'CDS'

        if feature:
            """Initialise class members(qualifiers) using SeqFeature object"""
            self._qualifiers = feature.qualifiers
            self.id = feature.id
            self.locus_tag = self._qualifiers.pop('locus_tag', [None])[0]
            self.product = self._qualifiers.pop('product', [None])[0]
            self.protein_id = self._qualifiers.pop('protein_id', [None])[0]
            self.gene = self._qualifiers.pop('gene', [None])[0]
            self.translation = self._qualifiers.pop('translation', [None])[0]
            self.notes = self._qualifiers.pop('note', [])
            self.EC_number = self._qualifiers.pop('EC_number', [])
            self.transl_table = self._qualifiers.pop('transl_table', [None])[0]
            self.source = self._qualifiers.pop('source', [None])[0]
            self.aSASF_note = self._qualifiers.pop('aSASF_note', [])
            self.aSASF_choice = self._qualifiers.pop('aSASF_choice', [])
            self.aSASF_scaffold = self._qualifiers.pop('aSASF_scaffold', [])
            self.aSASF_prediction = self._qualifiers.pop('aSASF_prediction', [])
            self.aSProdPred = self._qualifiers.pop('aSProdPred', [])
            self.db_xref = self._qualifiers.pop('db_xref', [])
            self.sec_met_predictions = self._qualifiers.pop('sec_met_predictions', [])
            self.location = feature.location
            if 'sec_met' in self._qualifiers:
                self.sec_met = self._map_sec_met_list_to_SecMetQualifier(self._qualifiers['sec_met'])
        else:
            self.location = f_location

    def get_cluster(self):
        """Returns the corresponding ClusterFeature"""
        return self.cluster

    def _map_sec_met_list_to_SecMetQualifier(self, sec_met_as_list):
        """Convert sec_met in list form to SecMetQualifier() form"""
        self._clustertype = None
        self._domains = None
        self._kind = None
        self._nrpspks = []
        self._asf_predictions = []
        for qualifier in sec_met_as_list:
            if qualifier.startswith('Type: '):
                self._clustertype = qualifier.split()[-1]
            elif qualifier.startswith('Kind: '):
                self._kind = qualifier.split()[-1]
            elif qualifier.startswith('Domains detected: '):
                qualifier = qualifier[18:]
                domains = qualifier.split(';')
                self._domains = []
                for domain in domains:
                    domain_name = domain.partition(" (")[0].replace(" ", "")
                    evalue = domain.partition("E-value: ")[2].partition(",")[0]
                    bitscore = domain.partition("bitscore: ")[2].partition(",")[0]
                    nr_seeds = domain.partition("seeds: ")[2].partition(")")[0]
                    sec_met_result = SecMetResult()
                    sec_met_result.query_id = domain_name
                    sec_met_result.evalue = evalue
                    sec_met_result.bitscore = bitscore
                    sec_met_result.nseeds = nr_seeds
                    self._domains.append(sec_met_result)
            elif qualifier.startswith('ASF-prediction: '):
                self._asf_predictions.append(qualifier)
            elif qualifier.startswith('NRPS/PKS '):
                self._nrpspks.append(qualifier)
        sec_met = SecMetQualifier(self._clustertype, self._domains, self._kind)
        if self._nrpspks:
            sec_met.nrpspks = self._nrpspks
        if self._asf_predictions:
            sec_met.asf_prediction = self._asf_predictions
        return sec_met

    def to_biopython(self):
        """Returns a Bio.SeqFeature.SeqFeature object with all its members"""
        new_CDS = SeqFeature(self.location, type=self.type, id=self.id)
        if not isinstance(self.sec_met, SecMetQualifier):
            raise ValueError('Invalid sec_met type')
        self._qualifiers['sec_met'] = self.sec_met.as_list()
        if self.locus_tag:
            self._qualifiers['locus_tag'] = [str(self.locus_tag)]
        if self.product:
            self._qualifiers['product'] = [str(self.product)]
        if self.protein_id:
            self._qualifiers['protein_id'] = [str(self.protein_id)]
        if self.gene:
            self._qualifiers['gene'] = [str(self.gene)]
        if self.translation:
            self._qualifiers['translation'] = [str(self.translation)]
        if self.notes:
            self._qualifiers['note'] = self.notes
        if self.EC_number:
            self._qualifiers['EC_number'] = self.EC_number
        if self.transl_table:
            self._qualifiers['transl_table'] = [str(self.transl_table)]
        if self.source:
            self._qualifiers['source'] = [str(self.source)]
        if self.db_xref:
            self._qualifiers['db_xref'] = self.db_xref
        if self.aSASF_choice:
            self._qualifiers['aSASF_choice'] = self.aSASF_choice
        if self.aSASF_note:
            self._qualifiers['aSASF_note'] = self.aSASF_note
        if self.aSASF_prediction:
            self._qualifiers['aSASF_prediction'] = self.aSASF_prediction
        if self.aSASF_scaffold:
            self._qualifiers['aSASF_scaffold'] = self.aSASF_scaffold
        if self.aSProdPred:
            self._qualifiers['aSProdPred'] = self.aSProdPred
        if self.sec_met_predictions:
            self._qualifiers['sec_met_predictions'] = self.sec_met_predictions
        new_CDS.qualifiers = self._qualifiers.copy()
        return [new_CDS]

    def __repr__(self):
        """A string representation of biopython CDS feature"""
        return repr(self.to_biopython()[0])


class SubCDSFeature(Feature):
    """A super class for CDS_motifFeature, PFAM_domain and aSDomain"""
    def __init__(self, f_location, feature):
        super(SubCDSFeature, self).__init__()
        self.asDomain_id = None
        self.aSTool = None
        self.detection = None
        self.database = None
        self.translation = None
        self.locus_tag = None
        self.label = None
        self.aSProdPred = []
        self.aSASF_choice = []
        self.aSASF_note = []
        self.aSASF_prediction = []
        self.aSASF_scaffold = []
        self._qualifiers = {}

        if feature:
            self._qualifiers = feature.qualifiers
            self.locus_tag = self._qualifiers.pop('locus_tag', [None])[0]
            self.translation = self._qualifiers.pop('translation', [None])[0]
            self.asDomain_id = self._qualifiers.pop('asDomain_id', [None])[0]
            self.aSTool = self._qualifiers.pop('aSTool', [None])[0]
            self.detection = self._qualifiers.pop('detection', [None])[0]
            self.database = self._qualifiers.pop('database', [None])[0]
            self.label = self._qualifiers.pop('label', [None])[0]
            self.aSASF_note = self._qualifiers.pop('aSASF_note', [])
            self.aSASF_choice = self._qualifiers.pop('aSASF_choice', [])
            self.aSASF_scaffold = self._qualifiers.pop('aSASF_scaffold', [])
            self.aSASF_prediction = self._qualifiers.pop('aSASF_prediction', [])
            self.aSProdPred = self._qualifiers.pop('aSProdPred', [])
            self.notes = self._qualifiers.pop('note', [])
            self.location = feature.location
            if 'score' in self._qualifiers:
                self.score = self._qualifiers['score'][0]
                del self._qualifiers['score']
            if 'evalue' in self._qualifiers:
                self.evalue = self._qualifiers['evalue'][0]
                del self._qualifiers['evalue']
        else:
            self.location = f_location
    #Check for a valid score qualifier before assigning
    def _get_score(self):
        try:
            return self.__score
        except:
            return None
    def _set_score(self, value):
        try:
            self.__score = float(value)
        except ValueError:
            raise ValueError('Invalid score value')
    score = property(_get_score, _set_score)

    #Check for a valid evalue qualifier before assigning
    def _get_evalue(self):
        try:
            return self.__evalue
        except:
            return None
    def _set_evalue(self, value):
        try:
            self.__evalue = float(value)
        except ValueError:
            raise ValueError('Invalid evalue value')
    evalue = property(_get_evalue, _set_evalue)

    def _get_feature_qualifiers(self):
        if self.locus_tag:
            self._qualifiers['locus_tag'] = [str(self.locus_tag)]
        if self.translation:
            self._qualifiers['translation'] = [str(self.translation)]
        if self.database:
            self._qualifiers['database'] = [str(self.database)]
        if self.evalue:
            self._qualifiers['evalue'] = [str(self.evalue)]
        if self.asDomain_id:
            self._qualifiers['asDomain_id'] = [str(self.asDomain_id)]
        if self.detection:
            self._qualifiers['detection'] = [str(self.detection)]
        if self.score:
            self._qualifiers['score'] = [str(self.score)]
        if self.label:
            self._qualifiers['label'] = [str(self.label)]
        if self.aSTool:
            self._qualifiers['aSTool'] = [str(self.aSTool)]
        if self.aSASF_choice:
            self._qualifiers['aSASF_choice'] = self.aSASF_choice
        if self.aSASF_note:
            self._qualifiers['aSASF_note'] = self.aSASF_note
        if self.aSASF_prediction:
            self._qualifiers['aSASF_prediction'] = self.aSASF_prediction
        if self.aSASF_scaffold:
            self._qualifiers['aSASF_scaffold'] = self.aSASF_scaffold
        if self.aSProdPred:
            self._qualifiers['aSProdPred'] = self.aSProdPred
        if self.notes:
            self._qualifiers['note'] = self.notes
        return self._qualifiers

class CDS_motifFeature(SubCDSFeature):
    """A CDS_motifFeature which subclasses Feature"""
    def __init__(self, f_location=None, feature=None):
        """Initialise a CDS_motifFeature
            param f_location: class 'Bio.SeqFeature.FeatureLocation/CompoundLocation'
            param feature: class 'Bio.SeqFeature.SeqFeature'
        """
        super(CDS_motifFeature, self).__init__(f_location, feature)
        self.motif = None
        self.type = 'CDS_motif'

        if feature:
            """Initialise class members(qualifiers) using SeqFeature object"""
            self.motif = feature.qualifiers.pop('motif', [None])[0]

    def to_biopython(self):
        """Returns a Bio.SeqFeature.SeqFeature object with all its members"""
        new_CDS_motif = SeqFeature(self.location, type=self.type)
        self._qualifiers = self._get_feature_qualifiers()
        if self.motif:
            self._qualifiers['motif'] = [str(self.motif)]
        new_CDS_motif.qualifiers = self._qualifiers.copy()
        return [new_CDS_motif]

    def __repr__(self):
        """A string representation of biopython CDS_motif feature"""
        return repr(self.to_biopython()[0])


class PFAM_domain(SubCDSFeature):
    """A PHAM_domain feature which subclasses Feature"""
    def __init__(self, f_location=None, feature=None):
        """Initialise a ClusterFeature
            param f_location: class 'Bio.SeqFeature.FeatureLocation/CompoundLocation'
            param feature: class 'Bio.SeqFeature.SeqFeature'
        """
        super(PFAM_domain, self).__init__(f_location, feature)
        self.domain = None
        self.description = None
        self.db_xref = []
        self.type = 'PFAM_domain'

        if feature:
            """Initialise class members(qualifiers) using SeqFeature object"""
            self.domain = feature.qualifiers.pop('domain', [None])[0]
            self.description = feature.qualifiers.pop('description', [None])[0]
            self.db_xref = feature.qualifiers.pop('db_xref', [])

    def to_biopython(self):
        """Returns a Bio.SeqFeature.SeqFeature object with all its members"""
        new_PFAM_domain = SeqFeature(self.location, type=self.type)
        self._qualifiers = self._get_feature_qualifiers()
        if self.domain:
            self._qualifiers['domain'] = [str(self.domain)]
        if self.description:
            self._qualifiers['description'] = [str(self.description)]
        if self.db_xref:
            self._qualifiers['db_xref'] = self.db_xref
        new_PFAM_domain.qualifiers = self._qualifiers.copy()
        return [new_PFAM_domain]

    def __repr__(self):
        """A string representation of biopython PFAM_domain feature"""
        return repr(self.to_biopython()[0])


class aSDomain(SubCDSFeature):
    """A aSDomain feature which subclasses Feature"""
    def __init__(self, f_location=None, feature=None):
        """Initialise a ClusterFeature
            param f_location: class 'Bio.SeqFeature.FeatureLocation/CompoundLocation'
            param feature: class 'Bio.SeqFeature.SeqFeature'
        """
        super(aSDomain, self).__init__(f_location, feature)
        self.domain = None
        self.domain_subtype = None
        self.specificity = []
        self.type = 'aSDomain'

        if feature:
            """Initialise class members(qualifiers) using SeqFeature object"""
            self.domain = feature.qualifiers.pop('domain', [None])[0]
            self.domain_subtype = feature.qualifiers.pop('domain_subtype', [None])[0]
            self.specificity = feature.qualifiers.pop('specificity', [])

    def to_biopython(self):
        """Returns a Bio.SeqFeature.SeqFeature object with all its members"""
        new_aSDomain = SeqFeature(self.location, type=self.type)
        self._qualifiers = self._get_feature_qualifiers()
        if self.domain:
            self._qualifiers['domain'] = [str(self.domain)]
        if self.domain_subtype:
            self._qualifiers['domain_subtype'] = [str(self.domain_subtype)]
        if self.specificity:
            self._qualifiers['specificity'] = self.specificity
        new_aSDomain.qualifiers = self._qualifiers.copy()
        return [new_aSDomain]

    def __repr__(self):
        """A string representation of the biopython aSDomain feature"""
        return repr(self.to_biopython()[0])


class ClusterFeature(Feature):
    """A ClusterFeature which subclasses Feature"""
    def __init__(self, f_location=None, feature=None):
        """Initialise a ClusterFeature
            param f_location: class 'Bio.SeqFeature.FeatureLocation/CompoundLocation'
            param feature: class 'Bio.SeqFeature.SeqFeature'
        """
        super(ClusterFeature, self).__init__()
        self.contig_edge = None
        self.detection = None
        self.products = []
        self._qualifiers = {}
        self.parent_record = None
        self.type = 'cluster'
        self.structure = None
        self.probability = None
        self.subclusterblast = None
        self.knownclusterblast = None
        self.clusterblast = None
        self.cdss = []

        if feature:
            """Initialise class members(qualifiers) using SeqFeature object"""
            self._qualifiers = feature.qualifiers
            self.contig_edge = self._qualifiers.pop('contig_edge', [None])[0]
            self.products = self._qualifiers.pop('product', [])
            self.structure = self._qualifiers.pop('structure', [None])[0]
            self.probability = self._qualifiers.pop('probability', [None])[0]
            self.subclusterblast = self._qualifiers.pop('subclusterblast', [])
            self.knownclusterblast = self._qualifiers.pop('knownclusterblast', [])
            self.clusterblast = self._qualifiers.pop('clusterblast', [])
            self.notes = self._qualifiers.pop('note', [])
            if self.notes:
                note_list_copy = self.notes[:]
                for  value in self.notes:
                    if value.startswith('Cluster number'):
                        note_list_copy.remove(value)
                    if value.startswith('Detection rule(s)'):
                        self.detection = value
                        note_list_copy.remove(value)
                self.notes = note_list_copy
            if 'cutoff' in self._qualifiers:
                self.cutoff = int(self._qualifiers['cutoff'][0])
            if 'extension' in self._qualifiers:
                self.extension = int(self._qualifiers['extension'][0])
            self.location = feature.location
        else:
            self.location = f_location

    #Check if cutoff is an integer before assigning
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

    #Check if extension is an integer before assigning
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
            return 0
        return self.parent_record.get_cluster_number(self)

    def get_CDSs(self):
        """Retruns a list of CDSFeature(s) which belong to this cluster"""
        return tuple(self.cdss)

    def to_biopython(self):
        """Returns a Bio.SeqFeature.SeqFeature object with all its members"""
        new_Cluster = SeqFeature(self.location, type=self.type)
        self._qualifiers['note'] = ["Cluster number: " + str(self.get_cluster_number())]
        if self.detection:
            self._qualifiers['note'].append(self.detection)
        if self.notes:
            self._qualifiers['note'].extend(self.notes)
        if self.cutoff:
            self._qualifiers['cutoff'] = [str(self.cutoff)]
        if self.extension:
            self._qualifiers['extension'] = [str(self.extension)]
        if self.products:
            self._qualifiers['product'] = self.products
        if self.contig_edge:
            self._qualifiers['contig_edge'] = [str(self.contig_edge)]
        if self.structure:
            self._qualifiers['structure'] = [str(self.structure)]
        if self.probability:
            self._qualifiers['probability'] = [str(self.probability)]
        if self.subclusterblast:
            self._qualifiers['subclusterblast'] = self.subclusterblast
        if self.knownclusterblast:
            self._qualifiers['knownclusterblast'] = self.knownclusterblast
        if self.clusterblast:
            self._qualifiers['clusterblast'] = self.clusterblast
        new_Cluster.qualifiers = self._qualifiers.copy()
        return [new_Cluster]

    def __repr__(self):
        """A string representation of biopython cluster feature"""
        return repr(self.to_biopython()[0])

class Record(object):
    """A record containing secondary metabolite clusters"""

    def __init__(self, seq_record=None):
        """Initialise a secondary metabolite record

        :param seq_record:  :class:`Bio.SeqRecord.SeqRecord` to read
        :type seq_record:   :class:`Bio.SeqRecord.SeqRecord`
        """
        self._record = seq_record
        self._modified_cds = []         #A list containing instances of CDSFeature
        self._modified_cluster = []     #A list containing instances of ClusterFeature
        self._modified_generic = []     #A list containing instances of GenericFeature
        self._modified_cds_motif = []   #A list containing instances of CDS_motifFeature
        self._modified_pfam_domain = [] #A list containing instances of PFAM_domain
        self._modified_asdomain = []    #A list containing instances of aSDomain
        self._cluster_number_dict = {}  #A dictionary to map clusters and their numbers

        if self._record:
            if not isinstance(self._record, SeqRecord):
                raise ValueError("SeqRecord should be an instance of 'Bio.SeqRecord.SeqRecord'")
            self.from_biopython(self._record)
        else:
            self._record = SeqRecord(Seq(""))

    @classmethod
    def from_file(cls, filename):
        """Initialise a record from a file

        :param string filename:    file name of the file to read
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
        if self._record:
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
        if self._record:
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
        if self._record:
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
        if self._record:
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
        if self._record:
            return self._record.annotations
        else:
            return {}
    def add_annotation(self, key, value):
        """Adding annotations in Record"""
        if not (isinstance(key, str) and (isinstance(value, str) or isinstance(value, list))):
            raise ValueError('Key and Value are not in right format')
        self._record.annotations[key] = value

    def __len__(self):
        """Return the length of the Biorecord"""
        return len(self._record)

    def get_clusters(self):
        """A list of secondary metabolite clusters present in the record"""
        return tuple(self._modified_cluster)
    def erase_clusters(self):
        """Erase all cluster features from the Record"""
        self._modified_cluster = []

    def get_CDSs(self):
        """A list of secondary metabolite CDS features present in the record"""
        return tuple(self._modified_cds)
    def erase_CDSs(self):
        """Erase all CDS features from the Record"""
        self._modified_cds = []

    def get_CDS_motifs(self):
        """A list of secondary metabolite CDS_motifs present in the record"""
        return tuple(self._modified_cds_motif)
    def erase_CDS_motifs(self):
        """Erase all CDS_motif features present in the Record"""
        self._modified_cds_motif = []

    def get_PFAM_domains(self):
        """A list of secondary metabolite PFAM_domains present in the record"""
        return tuple(self._modified_pfam_domain)
    def erase_PFAM_domains(self):
        """Erase all PFAM_domain features present in the Record"""
        self._modified_pfam_domain = []

    def get_aSDomains(self):
        """A list of secondary metabolite aSDomains present in the record"""
        return tuple(self._modified_asdomain)
    def erase_aSDomains(self):
        """Erase all aSDomain features present in the Record"""
        self._modified_asdomain = []

    def get_generics(self):
        """A list of secondary metabolite generics present in the record"""
        return tuple(self._modified_generic)
    def erase_generics(self):
        """Erase all generic features present in the Record"""
        self._modified_generic = []

    def to_biopython(self):
        """Returns a Bio.SeqRecord instance of the record"""
        new_record = self._record
        features = list(self.get_generics())
        features.extend(list(self.get_clusters()))
        features.extend(list(self.get_CDSs()))
        features.extend(list(self.get_CDS_motifs()))
        features.extend(list(self.get_aSDomains()))
        features.extend(list(self.get_PFAM_domains()))
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
        """Adds feature to appropriate lists"""
        if not isinstance(feature, Feature):
            raise TypeError("The argument is not an instance of 'Feature'")
        if isinstance(feature, ClusterFeature):
            clusters = self._modified_cluster
            index = find_new_cluster_pos(clusters, feature)
            clusters.insert(index, feature)
            feature.parent_record = self
            #Link cluster feature with its cds features
            self._update_cluster_cds_links(feature)
            for i, cluster in enumerate(clusters):
                self._cluster_number_dict[cluster] = i+1

        elif isinstance(feature, CDSFeature):
            self._modified_cds.append(feature)
            #Link cds feature with its cluster feature
            self._update_cluster_cds_links(feature)
        elif isinstance(feature, CDS_motifFeature):
            self._modified_cds_motif.append(feature)
        elif isinstance(feature, PFAM_domain):
            self._modified_pfam_domain.append(feature)
        elif isinstance(feature, aSDomain):
            self._modified_asdomain.append(feature)
        else:
            self._modified_generic.append(feature)

    def from_biopython(self, record):
        """Modifies _modified_features_* list with new Feature instances"""
        for feature in record.features:
            if feature.type == 'CDS':
                feature = CDSFeature(feature=feature)
                self._modified_cds.append(feature)
            elif feature.type == 'cluster':
                feature = ClusterFeature(feature=feature)
                feature.parent_record = self
                self._modified_cluster.append(feature)
                self._cluster_number_dict[feature] = self._modified_cluster.index(feature)+1
            elif feature.type == 'CDS_motif':
                feature = CDS_motifFeature(feature=feature)
                self._modified_cds_motif.append(feature)
            elif feature.type == 'PFAM_domain':
                feature = PFAM_domain(feature=feature)
                self._modified_pfam_domain.append(feature)
            elif feature.type == 'aSDomain':
                feature = aSDomain(feature=feature)
                self._modified_asdomain.append(feature)
            else:
                feature = GenericFeature(feature=feature)
                self._modified_generic.append(feature)
        return self

    def _update_cluster_cds_links(self, feature):
        """Link cluster and CDS features"""
        if isinstance(feature, ClusterFeature):
            clustercdsfeatures = []
            for cds in self.get_CDSs():
                if feature.location.start <= cds.location.start <= feature.location.end or \
                   feature.location.start <= cds.location.end <= feature.location.end:
                    clustercdsfeatures.append(cds)
                    cds.cluster = feature
            feature.cdss = clustercdsfeatures
        else:
            clusters = self.get_clusters()
            find_cluster_of_new_cds(clusters, feature)


class SecMetQualifier(list):
    """A SecMetQualifier class for sec_met qualifiers"""

    def __init__(self, clustertype=None, domains=None, kind=None):
        """Initialise a SecMetQualifier with the given attributes
            :param clustertype: an instance of str
            :param domains: a list of SecMetResult instance(s)
            :param kind: an instance of str
        """
        if clustertype and not isinstance(clustertype, str):
            raise TypeError('clustertype should be an instance of str')
        if domains and not isinstance(domains, list):
            raise TypeError('domains should be an instance of list')
        if kind and not isinstance(kind, str):
            raise TypeError('kind should be an instance of str')
        self.clustertype = clustertype
        self.domains = domains
        self.kind = kind
        self.nrpspks = []
        self.asf_predictions = []
        super(SecMetQualifier, self).__init__()

    def __len__(self):
        """Return length of the secmet qualifier"""
        count = 0
        if self.clustertype:
            count += 1
        if self.domains:
            count += 1
        if self.kind:
            count += 1
        if self.nrpspks:
            count += len(self.nrpspks)
        if self.asf_predictions:
            count += len(self.asf_predictions)
        return count

    def __repr__(self):
        """A string representation of the list of sec_met qualifier"""
        return str(self.as_list())

    def __iter__(self):
        if self.clustertype:
            yield "Type: %s" % self.clustertype
        if self.domains:
            yield "Domains detected: " + "; ".join(map(str, self.domains))
        if self.kind:
            yield "Kind: %s" % self.kind
        if self.nrpspks:
            for nrps in self.nrpspks:
                yield nrps
        if self.asf_predictions:
            for asf in self.asf_predictions:
                yield asf

    def as_list(self):
        """Returns sec_met qualifier in a list"""
        self._sec_met = []
        for qual in self:
            self._sec_met.append(qual)
        return self._sec_met

class SecMetResult():
    def __init__(self, res=None, nseeds=None):
        self.query_id = None
        self.evalue = None
        self.bitscore = None
        self.nseeds = None
        if res and nseeds:
            self.query_id = res.query_id
            self.evalue = res.evalue
            self.bitscore = res.bitscore
            self.nseeds = nseeds

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "{} (E-value: {}, bitscore: {}, seeds: {})".format(
            self.query_id, self.evalue, self.bitscore, self.nseeds)
