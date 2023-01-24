import time
import pickle
from typing import List

from crunch_classes import CandidateWithOffTargets


def remove_duplicates(list_of_candidates: List[CandidateWithOffTargets]):
    """
    This function removes duplicated sgRNAs with identical sequences.
    :param list_of_candidates: A list of CandidateWithOffTargets objects
    :return: A list of all CandidateWithOffTargets objects from all subgroups, without duplicates.
    """
    sequence_to_candidate_dict = {}
    for candidate in list_of_candidates:
        if candidate.seq not in sequence_to_candidate_dict:
            sequence_to_candidate_dict[candidate.seq] = candidate
            continue
        sequence_to_candidate_dict[candidate.seq] = choose_candidate_from_duplicates(
            sequence_to_candidate_dict[candidate.seq], candidate)
    return [sequence_to_candidate_dict[key] for key in sequence_to_candidate_dict]


def choose_candidate_from_duplicates(candidate_1: CandidateWithOffTargets, candidate_2: CandidateWithOffTargets):
    """
    This function picks a candidate between two duplicated sgRNA candidates.
    :param candidate_1: a CandidateWithOffTargets object.
    :param candidate_2:a CandidateWithOffTargets object.
    :return: The candidate with the higher cut expectation. If both candidates have the same scores,
    the candidate that belongs to the smaller (descendant) subgroup will be returned.
    """
    assert candidate_2.seq == candidate_1.seq
    # Pick the candidate with the higher cut expectation.
    if candidate_1.cut_expectation > candidate_2.cut_expectation:
        return candidate_1
    elif candidate_1.cut_expectation == candidate_2.cut_expectation:
        # Pick the candidate from the smaller (descendant) subgroup.
        if len(candidate_1.subgroup.genes_lst) < len(candidate_2.subgroup.genes_lst):
            return candidate_1
        return candidate_2
    return candidate_2


def remove_overlapping(list_of_candidates: List[CandidateWithOffTargets], overlap_threshold: int = 10):
    """
    This function removes sgRNA candidates that target the genome in overlapping regions. It always removes the
    worse-scored candidates.
    :param list_of_candidates: A list of CandidateWithOfftarget objects.
    :param overlap_threshold:
    :return: A list of CandidateWithOfftarget objects without overlaps.
    """
    genes_covered_to_candidates_dict = {}
    for candidate in list_of_candidates:
        genes_covered = tuple(sorted(candidate.genes_score_dict.keys()))
        # If the candidate covers a new set of targeted genes, it doesn't overlap with any of the targets already in the
        # dictionary.
        if genes_covered not in genes_covered_to_candidates_dict:
            genes_covered_to_candidates_dict[genes_covered] = [candidate]
            continue
        # Else, there are other candidates that target the same genes, go over them and check for overlaps.
        overlapping = False
        for i, other_candidate in enumerate(genes_covered_to_candidates_dict[genes_covered]):
            overlapping = is_overlapping(candidate, other_candidate, overlap_threshold)
            if overlapping:
                # If the cut expectation of the new candidate is higher than the old overlapping candidate, replace the
                # old candidate with the new.
                if candidate.cut_expectation > other_candidate.cut_expectation:
                    genes_covered_to_candidates_dict[genes_covered][i] = candidate
                    break
        # If no overlap was found, add the candidate to the dictionary.
        if not overlapping:
            genes_covered_to_candidates_dict[genes_covered].append(candidate)
    output_list = []
    for val in genes_covered_to_candidates_dict.values():
        output_list += val
    return output_list


def create_sorted_targets_list(list_of_targets):
    """
    This function takes a list of targets of a specific gene, as given by each value in CandidateWithOffTargets.targets_dict[gene],
    and sorts them according to the strand and coordinates.
    :param list_of_targets: a nested list of lists, where each element represents a genomic target.
    :return: a nested list sorted based on strand and coordinates
    """
    return sorted(list_of_targets, key=lambda x: (x[4], x[3]))


def is_overlapping(candidate_1: CandidateWithOffTargets, candidate_2: CandidateWithOffTargets,
                   overlap_threshold: int = 10) -> bool:
    """
    Checks whether two sgRNA candidates overlap in all targets.
    :param candidate_1: a CandidateWithOfftarget object
    :param candidate_2: a CandidateWithOfftarget object
    :param overlap_threshold: a threshold for the minimum number of overlapping positions.
    :return: If the candidates overlap at all targets, return True. Otherwise, return False.
    """
    for gene in candidate_1.targets_dict:
        targets_list_1 = create_sorted_targets_list(candidate_1.targets_dict[gene])
        targets_list_2 = create_sorted_targets_list(candidate_2.targets_dict[gene])
        if len(targets_list_1) != len(targets_list_2):
            # If candidate 1 has different number of targets than candidate 2, they do not overlap.
            return False
        for target_1, target_2 in zip(targets_list_1, targets_list_2):
            range_target_1, range_target_2 = create_target_range(target_1), create_target_range(target_2)
            if range_target_1 == range_target_2:
                continue
            overlap = range(max(range_target_1.start, range_target_2.start),
                            min(range_target_1.stop, range_target_2.stop))
            if len(overlap) >= overlap_threshold:
                continue
            # If no overlaps were found between the targets in the gene, the candidates do not fully overlap
            return False
    return True


def create_target_range(target):
    """
    This function finds the coordinate range of a given target.
    :param target: An OffTarget object
    :return: The coordinate range for the target.
    """
    return range(target[3], target[3] + 20)
