import pickle
from typing import List, Dict

import globals
from crunch_classes import CandidateWithOffTargets, OffTarget


def filter_main(list_of_candidates: List[CandidateWithOffTargets], path_to_pickle: str, feature_score_dict: Dict,
                ignore_regions: set, output_name ="crispys_output_crunch", genes_to_ignore_set = set()):
    """
    wrapper for filter candidates base on off-target score
    Args:
        genes_to_ignore_set: a set of genes to ignore when looking at off targets.
        crispys_output_name: The output name for the files created by the filter.
        path_to_pickle: path to a folder where the results of off-target tagging of crispys output
        feature_score_dict: a dictionary that give the region name and the score cutoff, {'region': 'score'}
        ignore_regions: a set of genomic region that will be excluded from filtering. for example a pseudogene is a region we dont want to filter if it is hit by sgRNA.
        list_of_candidates: A list of CandidateWithOfftarget objects.

    Returns: the function will write to the result folder, 2 pickle files.
            one is a list of candidate that passed the filtering and the second is of candidate that do not pass.
            it also return the list of remain candidates

    """

    removed_candidates, remaining_candidate = filter_candidates_list(list_of_candidates, feature_score_dict,
                                                                  ignore_regions, genes_to_ignore_set)
    print(f"Feature scores used for off-target evaluation: {globals.feature_score_dict}")
    # write removed candidates to pickle
    with open(f"{path_to_pickle}/{output_name}_removed_candidates.p", "wb") as f:
        pickle.dump(removed_candidates, f)

    return remaining_candidate


def filter_candidates_list(list_of_candidates: List[CandidateWithOffTargets], feature_score_dict: dict,
                           ignore_regions: set = set(), genes_to_ignore_set: set = set()) -> tuple:
    """
    This function take a list of CandidateWithOffTargets objects and for each one go over all off-targets and filter the candidate based on
    the genomic region and score cutoff in the dictionary, any region that apper in the key will be filtered based on the score cutoff
     (candidate with off-target that has a regions with score above the cutoff will be discarded)
     It can also take a list of regions names that needs to be ignored from filtering.
    Args:
        genes_to_ignore_set: a set of genes to ignore when looking at off targets.
        list_of_candidates: list of CandidateWithOffTargets objects
        feature_score_dict: a dictionary that give the region name and the score cutoff, {'region': 'score'}
        ignore_regions: a set of genomic region that will be excluded from filtering.

    Returns:
        A tuple of two lists. first: a list of candidates that are discarded.
                              second: list of candidate that are kept.

    """
    removed_candidates = []
    remaining_candidates = []
    for candidate in list_of_candidates:
        keep_candidate = True
        for off_target in candidate.off_targets_list:
            # If the OffTarget cuts in one of the family genes, move to the next potential off-target.
            if not check_if_off_target(candidate, off_target, genes_to_ignore_set):
                continue
            check_if_passes_filter(off_target, ignore_regions, feature_score_dict)
            if not off_target.pass_filter:
                keep_candidate = False
                # add to removed list (if not already there)
                removed_seq = [can.seq for can in removed_candidates if removed_candidates]
                if candidate.seq not in removed_seq:
                    removed_candidates.append(candidate)
                # alert the user that there is a gene targeted which we couldn't get its name.
                if len(off_target.genes_covered) == 0:
                    print(f"Off target in gene without annotation: {off_target}\n")
                #IMPORTANT. now the code will exit in order to not screen more off-targets after the candidate is
                # rejected, that means that some off target that do not pass the threshold will still be marked with
                # off_target.pass_filter = True, to avoid that uncomment the 'break' statement
                break
        if keep_candidate:
            remaining_candidates.append(candidate)
        # if none of the off-targets give undesired results, keep the sgRNA.

    return removed_candidates, remaining_candidates


def check_if_off_target(candidate: CandidateWithOffTargets, off_target: OffTarget, genes_to_ignore_set: set = set()) -> bool:
    """
    This function checks if the off-target is located on a gene that should be targeted by the sgRNA candidate.
    genes_to_ignore_set: a set of genes to ignore when looking at off targets.
    :return: True if the target is located in a targeted gene, False otherwise. Update family_genes_to_off_targets_dict
    for the candidate if needed.
    """
    # check if the off-target hit a gene in candidate.genes_score_dict.
    family_genes_set = set(candidate.genes_score_dict.keys()).union(genes_to_ignore_set)
    if any(gene_name in family_genes_set for gene_name in off_target.genes_covered):
        off_target.in_family_gene = True
        # add the gene_name:off-target to the candidate.family_genes dictionary.
        gene_name = "".join(off_target.genes_covered)
        if gene_name not in candidate.family_genes_to_off_targets_dict.keys():
            candidate.family_genes_to_off_targets_dict[gene_name] = [off_target]
        else:
            candidate.family_genes_to_off_targets_dict[gene_name].append(off_target)
        return False
    return True


def check_if_passes_filter(off_target: OffTarget, ignore_regions: set, feature_score_dict: Dict):
    """
    This function checks if the off target passes the thresholds in feature_score_dict.
    :param removed_candidates:
    :param feature_score_dict:
    :param ignore_regions:
    :param candidate:
    :param off_target:
    :return:
    """
    ""
    # check if the off target score is higher than the threshold for any (all) regions
    if "Any" in feature_score_dict:
        if float(off_target.score) > float(feature_score_dict["Any"]):
            off_target.pass_filter = False
    # check if there is a feature that needs to be ignored, if so skip to next off-target.
    if any(region in off_target.genomic_regions for region in ignore_regions):
        return
    # go over the genomic regions covered by the off-target
    for region in off_target.genomic_regions:
        # if the off-target score is above the threshold specified in the dictionary, filter out the candidate.
        if region in feature_score_dict and (float(off_target.score) > float(feature_score_dict[region])):
            off_target.pass_filter = False




if __name__ == '__main__':

    feature_to_score_dict = {"gene": 0, "exon": 0}

    regions_to_ignore = {"pseudogene"}

    filter_main("/groups/itay_mayrose/udiland/crispys_off_target/crispys_out_HOM04D000944", feature_to_score_dict,
                regions_to_ignore)

    # filter_main("/groups/itay_mayrose/udiland/crispys_off_target/crispys_out_HOM04D000944", feature_to_score_dict)

    with open(
            "/groups/itay_mayrose/udiland/crispys_off_target/crispys_out_HOM04D000944/res_with_off_targets_filtered.p",
            "rb") as f:
        pass_filter = pickle.load(f)

    for sgrna_candidate in pass_filter[0].candidates_list:
        for off in sgrna_candidate.off_targets_list:
            print(off.genes_covered)
            print(off.genomic_regions)
