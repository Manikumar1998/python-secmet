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


    @property
    def id(self):
        """Pass through to seq_record object if available"""
        if self._record is not None:
            return self._record.id
        else:
            return "NO_ID_ASSIGNED"

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
