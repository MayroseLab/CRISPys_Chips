from Candidate import Candidate
from typing import List, Dict

class CandidateWithOffTargets(Candidate):
    """
    This class inherit the Candidate class and adds some off-target attributes
    """
    def __init__(self, candidate):
        super().__init__(candidate.seq, candidate.cut_expectation, candidate.genes_score_dict, candidate.targets_dict)
        self.off_targets_list = []
        self.family_genes_to_off_targets_dict = dict()
        self.subgroup = None


    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return self.__str__()


class OffTarget:
    """
    This class contains off-target information for an sgRNA candidate, it is intend to be part of a
    CandidateWithOffTarget object as an item of the off_target_list
    """
    def __init__(self, seq: str, chromosome: str, start_position: int, strand: str, number_of_mismatches: int):
        self.seq = seq
        self.chromosome = chromosome
        self.start_position = start_position
        self.strand = strand
        self.number_of_mismatches = number_of_mismatches
        self.genes_covered = set()
        self.genomic_regions = set()
        self.score = -1
        self.pass_filter = True
        self.in_family_gene = False

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return f"[{self.seq}, {str(self.chromosome)}, {self.start_position}, {self.strand}, {self.number_of_mismatches}, {round(float(self.score), 4)}, {self.genomic_regions}, {self.genes_covered}]"

