import os
import time
import pickle
import copy
from typing import List, Dict
from Bio.Seq import reverse_complement
import random
from off_target_search import off_target_search_crispritz_main
from filter_off_targets import filter_main
import argparse
import globals
from make_tree_display_CSV import create_output_multiplex
import Candidate
from SubgroupRes import SubgroupRes
from crunch_classes import CandidateWithOffTargets
random.seed(42)
"""
This module will takes CRISPys output (list of SubgroupRes) and will output group of n guides that will target the 
maximum number of the genes in a group of genes (internal node).
It will first select the 'best' guide, that is, the one that will target the most genes with maximum expectation score, 
than, it will select another guide that will complement the target to capture more genes that are not capture with the 
previous guide, this will continue until we get the number of guide we specify (usually the multiplex is of 2 guides).
Next, this process will repeat but without selecting the 'additional' guides, meaning, it will keep the first ('best') 
guide and remove from the list of all guides that completed the multiplex in previous step, and will find new ones to add to the 'bset' guide. 
the 'best' one, this argument is controlled with the 'n_multiplx' parameter.
each multiplex (list of guides for one vector) is stored in a SubgroupRes object and a group of them are stored in a BestSgGroup object.
A representation of the BestSgGroup:


                        --------------
                        | multiplx 1 | each multiplex is a SubgroupRes object with n candidates (specify in the n_sgrnas parameter) 
                        | multiplx 2 |
BestSgGroup object ->   |    .       |
                        |    .       |
                        | multiplx n |
                        --------------
each BestSgGroup has 'best candidate' sttribute that show the details of the guide that is common to all multiplex.

for each CRISPys output of gRNAs that target group of genes we can have multiple BestSgGroup objects each one with its 
own 'best' guide, so the next step is to remove the 'best guide' we selected from the CRISPys results and repeat the
first part again to create the next group of multiplex, this parameter is defined in 'number_of_groups'
so we get:

   BestSgGroup object 1    BestSgGroup object 2   .   .   .  BestSgGroup object n
   --------------          --------------                     --------------
   | multiplx 1 |          | multiplx 1 |                     | multiplx 1 |
   | multiplx 2 |          | multiplx 2 |                     | multiplx 2 |
   |    .       |          |    .       |        .    .   .   |    .       |
   |    .       |          |    .       |                     |    .       |
   | multiplx n |          | multiplx n |                     | multiplx n |
   --------------          --------------                     --------------
The described above is applied to each subgroup (represent internal node in the gene tree) of CRISPys results and
the output object is a dictionary that each key is the internal node name and the value is a dictionary that each key is
the sequence of the 'best' gRNA and the value is the BestSgGroup object of that 'best guide' 
(in each object are all multiplexes for that 'best guide') 

"""


class BestSgGroup:
    """
    This class is made to store a collection of subgroups object, such that each contains group of candidates that where
    selected fo multiplexing.
    The multiplexing group is designed so that each subgroup contain one candidate that is the same
     (referred to as 'best candidate') and another (one or more) that are different in each group this why the class is
     called BestSgGroup and the 'best_candidate' attribute will store the sequence of the best candidate.
    """

    def __init__(self, best_candidate: Candidate.Candidate = None, subgroups_list: List = None,
                 all_candidates=None):
        if all_candidates is None:
            all_candidates = list()
        self.best_candidate = best_candidate
        self.subgroups_list = subgroups_list
        self.all_candidates = all_candidates

    def __str__(self):
        return f"{self.best_candidate} , {self.subgroups_list}"


def check_if_one_gene(subgroup_list: List):
    """
    This function check if crispys result contain output for only one gene
    and if so will quit the program
    Args:
        subgroup_list: list of SubgroupRes

    Returns:
    """
    # if the results are only singletons check how many genes in the results and if only one gene targeted exit the program
    number_of_genes_covered = len(set([str(candidate.genes_score_dict.keys()) for subgroup in subgroup_list for candidate in subgroup.candidates_list]))
    if number_of_genes_covered == 1:
        print(f"Only one gene targeted in {subgroup_list[0].family_name}")
        exit()



def add_singletons_to_subgroup(subgroup_list: list, number_of_singletons: int = 5) -> list:
    """
    This function takes the output of crispys that contain subgroupRes object with singletons (targeting one gene)
    and add the candidates of a singletons results to the object holding candidate for multiple genes (internal node in the gene tree)
    which one of the is the gene targeted by the singleton candidates.
    notice that the function add singltones by looking at the genes targeted in the internal node and not all the gene in the internal node,
    if you want to add singletons targeting all genes in the internal node change 'subgroup.genes_lst' to 'subgroup.genes_in_node'
    Args:
        subgroup_list: list of subgroupRes objects (output of crispys)
        number_of_singltones: the number of singletons from each gene to add to the subgroups of internal node

    Returns:
            A list of subgroupRes objects each one contain candidates for internal node in the gene tree (with singletons)
    """
    singletons_dict = {}
    subgroup_list_no_singles = []
    # create new list with only subgroups targetiong multiple genes and store the singletons in a dictionary of genes_name:candidates_list
    for subgroup in subgroup_list:
        if len(subgroup.genes_in_node) == 1:
            singletons_dict[subgroup.genes_lst[0]] = subgroup.candidates_list
        else:
            subgroup_list_no_singles.append(subgroup)

    # sort singletons by on target score
    singletons_dict = {k: sorted(v, key=lambda x: x.on_target_score, reverse=True) for k, v in singletons_dict.items()}

    # go over the list of subgroups with no singletons and add the singleton candidate if it is not exists
    # the number of singletons to add for each gene is defined in the 'number_of_singltones' variable
    for subgroup in subgroup_list_no_singles:
        singletons2add = []
        # make a list of candidates sequences from the no-singletons subgroup
        subgroup_cand_seqs = [can.seq for can in subgroup.candidates_list]
        # go over the singletons dictionary and check if the singleton target a gene of the subgroup
        for single_gene in singletons_dict.keys():
            if single_gene in subgroup.genes_in_node:
                # for 'number_of_singltones' times go over each candidate in the list and if it is not in the subgroup,
                # add it to the list of singltones
                i = 0
                for single_candidate in singletons_dict[single_gene][0: number_of_singletons]:
                        if i == number_of_singletons:
                            break
                        if single_candidate.seq not in subgroup_cand_seqs:
                            # add to the selected singleton the subgroup attribute of the subgroup he will be added to
                            single_candidate2add = copy.deepcopy(single_candidate)
                            single_candidate2add.subgroup = subgroup
                            # accumulate singletons to later add to subgroup
                            singletons2add.append(single_candidate2add)
                            i += 1
        # shuffle the singletons list and add it to subgroup candidate list
        random.shuffle(singletons2add)
        subgroup.candidates_list.extend(singletons2add)
    return subgroup_list_no_singles


def create_candidatewithofftargets_object(candidate, subgroup):
    """
    This function takes a Candidates object and returns a CandidateWithOffTargets object.
    :param candidate: a Candidate object.
    :param subgroup: a SubgroupRes object.
    :return: a CandidateWithOffTargets object
    """
    new_candidate = CandidateWithOffTargets(candidate)
    new_candidate.subgroup = SubgroupRes(sorted(subgroup.genes_lst), [], subgroup.name, subgroup.genes_in_node,
                                         family_name=subgroup.family_name)
    new_candidate.on_target_score = candidate.on_target_score
    return new_candidate

def concat_candidates_from_subgroups(list_of_subgroup_results: List[SubgroupRes]) -> List:
    """
    This function unifies all subgroup results into a single list of CandidateWithOffTarget objects.
    :param list_of_subgroup_results: list of Subgroup.
    :return: a list of CandidateWithOffTargets objects.
    """
    list_of_candidates = []
    for subgroup in list_of_subgroup_results:
        if subgroup.candidates_list:
            for candidate in subgroup.candidates_list:
                new_candidate = create_candidatewithofftargets_object(candidate, subgroup)
                list_of_candidates.append(new_candidate)
    return list_of_candidates

def remove_candidates_with_restriction_site(list_of_candidates: List[CandidateWithOffTargets],
                                            restriction_site: str = "None"):
    """
    This function returns a list of candidates which don't contain a given restriction site.

    Returns:

    """
    if restriction_site == "None":
        return list_of_candidates
    return [candidate for candidate in list_of_candidates if restriction_site not in candidate.seq and
            restriction_site not in reverse_complement(candidate.seq)]

def recreate_subgroup_lst(list_of_candidates_filtered: List, list_of_subgroup_results) -> List:
    """
    This function take a list of candidates (that passed off-targets filtering) and output a list os subgroupres objects
    each contain the candidate of an internal node of the gene tree (the same structure of crispys output)
    Args:
        list_of_candidates_filtered: list of candidate objects

    Returns:
        list of SubgroupRes
    """

    subgroups_dict = {}
    for candidate in list_of_candidates_filtered:
        if tuple(candidate.subgroup.genes_in_node) not in subgroups_dict.keys():
            subgroups_dict[tuple(candidate.subgroup.genes_in_node)] = [candidate]
        else:
            subgroups_dict[tuple(candidate.subgroup.genes_in_node)].append(candidate)
    # make a list of subgroupRes objects
    subgroup_lst = [SubgroupRes(list(genes), can_lst, can_lst[0].subgroup.name, list(genes)) for genes, can_lst in subgroups_dict.items()]
    # add subgroup that has no candidates (will be added with empty candidates_list)
    for subgroup in list_of_subgroup_results:
        if not subgroup.candidates_list:
            subgroup_lst += [SubgroupRes(subgroup.genes_lst, [], subgroup.name, subgroup.genes_in_node)]

    return subgroup_lst


def subgroup2dict(subgroup: SubgroupRes) -> Dict:
    """
    This function take a subgroup object and output a dictionary of sequence:candidate
    Args:
        subgroup: subgroups objects
    Returns: a dictionary of candidates
    """
    candidates_dict = {}
    for candidate in subgroup.candidates_list:
        # if the candidate is in the dictionary replace it with the same one that have higher cut expectation (if exist)
        if candidate.seq not in candidates_dict or candidate.cut_expectation > candidates_dict[candidate.seq].cut_expectation:
            candidates_dict[candidate.seq] = candidate
    return candidates_dict


def select_secondary(primery_grna: CandidateWithOffTargets, candidates_list: List, lambda_coef: float = 0.9):
    """
    This function get the primary gRNA and a list of other gRNAs and return the gRNA that is the best match for the
    primary to target the genes in the internal node
    Args:
        primery_grna: gRNA that has been selected to be the primary in a multiplex
        candidates_list: list of gRNAS
        lambda_coef: a coefficient to multiply by the score of a gRNA, this is used in the probability calculation

    Returns:
        the gRNA that best complement the primary gRNA to target all genes
    """
    # get all genes in node
    genes = primery_grna.subgroup.genes_in_node
    max_cut_probability = 0
    chosen_grna = None
    # go over every gRNA and calculate its probability to cut all the gene of in node together with the primary gRNA
    # and return the secondary gRNA that has the highest probability
    for grna in candidates_list:
        current_grna_prob = 0
        for gene in genes:
            current_grna_prob += get_gene_score(primery_grna, grna, gene, lambda_coef)
        if current_grna_prob > max_cut_probability:
            max_cut_probability = current_grna_prob
            chosen_grna = grna
    if max_cut_probability == 0:
        return
    return chosen_grna


def get_gene_score(primary_grna: CandidateWithOffTargets, grna: CandidateWithOffTargets, gene: str,
                   lambda_coef: float = 0.9) -> float:
    """
    This function calculate the probability of two gRNAs to cut a gene
    specifically its calculate the probability of NOT targeting the gene with both guides which is the multiplication of
    the (1 - score) for guide 1 with (1 - score) of  guide 2.
    Than it returns 1 - (p of not targeting both) which is the probability to target at least one of them.
    Args:
        primary_grna: gRNA1 - the guide that was selected as primary
        grna: gRNA2 - the guide that is tested to be matched with the primary
        gene: the name of the gene to use to calculate the probability
        lambda_coef: a scaling factor of the score for each target it is used

    Returns:
            the probability of the two guides to target the gene.
    """
    prob = 1 - ((1 - (lambda_coef * primary_grna.genes_score_dict.get(gene, 0)))*(1 - (lambda_coef * grna.genes_score_dict.get(gene, 0))))
    return prob


def check_overlap_positions(candidate: Candidate, can_pos_dict=None):
    """
    This function check if candidates overlap, it used in two stages: 1) check for overlap candidates inside multiplex
    and 2) check if 'best candidate' overlap with previous 'best' that been chosen
    Args:
        candidate: A candidate object
        can_pos_dict: positions dictionary of gene:positions (to compare with)
        check_bestgroup_overlap: a flag to know if comapring 'best' candidates
    Returns:

    """
    # if no position dictionary supplied create and return such dictionary
    if not can_pos_dict:
        can_pos_dict = {candidate.seq: dict()}
        for gene, target in candidate.targets_dict.items():
            # take only the positions in the best match of the gene (the first item in the list)
            pos_set = {(target[0][3], target[0][4])}
            can_pos_dict[candidate.seq][gene] = pos_set
        return can_pos_dict

    if can_pos_dict:
        for can_seq in can_pos_dict.keys():
            overlaps = []
            for gene in candidate.targets_dict:
                # get candidate positions
                target = candidate.targets_dict[gene]
                can_pos_set = {(target[0][3], target[0][4])}
                try:
                    if can_pos_dict[can_seq][gene] == can_pos_set:
                        overlaps.append(True)
                    else:
                        overlaps.append(False)
                # if the gene is not in the positional dictionary they are not fully overlap
                except KeyError:
                    overlaps.append(False)

            if all(overlaps):
                return True
        # add the new target to dictionary
        can_pos_dict[candidate.seq] = dict()
        for gene in candidate.targets_dict.keys():
            target = candidate.targets_dict[gene][0]
            pos_set = {(target[3], target[4])}
            can_pos_dict[candidate.seq][gene] = pos_set

        return False

def check_duplicates(selected_cand_dict: dict, guide_seq_dict: Dict=dict()):
    """
    This function is used to check for duplicate multiplexes
    (it uses the global variable 'guide_seq_dict' to store the sequence and gene count of each multiplex)
    Args:
        selected_cand_dict: a dictionary of a multiplex with seq:candidate

    Returns:
        for the first call of the function it returns False while creating a dictionary of mulitiplex seqs tuple and
        the number of genes targeted (seq1,seq2,..): n_genes

        for the other calls will return True if the multiplex should be skipped because its duplicate with higehr
        number of genes (higher internal node) or False if it did not (keep the multiplx)
    """
    if not guide_seq_dict:
        seqs = tuple(sorted([seq for seq in selected_cand_dict.keys()]))
        # count the number of genes target in the multiplex (include repetition)
        n_genes_targeted = sum([len(guide.subgroup.genes_in_node) for guide in selected_cand_dict.values()])
        # create a dictionary of guides seq and number of genes targeted
        guide_seq_dict[seqs] = n_genes_targeted
        return False
    # get the sequences of the multiplex and check if they are in the dictionary of previously selected multiplex
    multiplx_seqs = {tuple(sorted([seq for seq in selected_cand_dict.keys()]))}
    existing_seqs = {seqs for seqs in guide_seq_dict.keys()}
    # if the same multiplx has been selected count the number of genes in node for each
    if multiplx_seqs.intersection(existing_seqs):
        #count the genes it node of current guide
        sum_genes_targeted_new = sum([len(guide.subgroup.genes_in_node) for guide in selected_cand_dict.values()])
        # get the previous gene count
        sum_genes_targeted_multi_exist = guide_seq_dict[tuple(multiplx_seqs)[0]]
        if sum_genes_targeted_multi_exist < sum_genes_targeted_new:
            return True
        else:
            return False
    # if the multiplex hasn't been called add its sequences and gene count to the dictionary
    guide_seq_dict[tuple(multiplx_seqs)[0]] = sum([len(guide.targets_dict.keys()) for guide in selected_cand_dict.values()])
    return False


def choose_multiplx(subgroup: SubgroupRes, n_sgrnas: int = 2, best_candidate: Candidate = None,
                    pos_dict: Dict = None):
    """
    This function takes SubgroupRes object and returns n guides that will target
     as many genes as possible. this is the the function that produce the multiplex with n_sgrnas guides targeting the
     most genes in the subgroup
     when it is run for the first time on subgroup it is run without the a 'best_candidate' argument and it retruns
     the first multiplex with the 'Best sgRNA' as subgroup object (an object used in CRISPys to store results of
     internal node, used here for multiplex).
     When it is run subsequently with the 'Best_candidate' argument it finds the rest of the gRNA for the multiplex and
     retruns the multiplex (as subgroup object)
    Args:
        subgroup: a subgroup obhect conatining a list of candidates
        n_sgrnas: number of guide to output
        best_candidate: a 'best' candidate that was already chosen
        pos_dict:

    Returns: subgroup object containing a list of candidates, Candidate object containing the 'best candidate'

    """

    # initate a dict that will store the gRNAs selected
    selected_candidates = {}

    # select best guide according to the cut_expectation score
    if not best_candidate:
        if len(subgroup.candidates_list) == 1:
            return None
        best_candidate = sorted(subgroup.candidates_list, key=lambda x: x.cut_expectation, reverse=True)[0]
        # calculate the guide position
        pos_dict = check_overlap_positions(best_candidate)
        # remove the best candidate from the subgroup
        subgroup.candidates_list.remove(best_candidate)
    # store the selected 'best' guide in a dictionary
    selected_candidates[best_candidate.seq] = best_candidate

    # select the rest (other than the best) of the candidates for the amount specified in n_sgrnas
    i = 1
    while i < n_sgrnas:
        # check if no candidates left
        if not subgroup.candidates_list:
            return None
        # select the secondary gRNAs
        candidate = select_secondary(best_candidate, subgroup.candidates_list)
        # calculate the guide position
        skip_candidate = check_overlap_positions(candidate, pos_dict)

        if skip_candidate:
            subgroup.candidates_list.remove(candidate)
            continue

        # store the selected guide in a dictionary
        selected_candidates[candidate.seq] = candidate

        # catch the end of the iteration to check if the multiplex is duplicate
        if (i == (n_sgrnas - 1)):
            if check_duplicates(selected_candidates):
                # del (candidates_dict[candidate.seq])
                del (selected_candidates[candidate.seq])
                continue

        i += 1
    # make output to a subgroup list
    cand_list = [can for can in selected_candidates.values()]
    genes = []
    for can in selected_candidates.values():
        genes += can.genes_score_dict.keys()
    genes_lst = genes
    # make a tuple of candidates sequence and use it as a name for the subgroup
    name = tuple(cand.seq for cand in cand_list)
    subgroup = SubgroupRes(list(set(genes_lst)), cand_list, name, subgroup.genes_in_node)
    return subgroup, best_candidate, pos_dict


def get_best_groups(subgroup: SubgroupRes, n_multiplx: int, n_sgrnas: int = 2):
    """
    This function takes crispys output as SubgroupRes object and finds n_multiplx of n_sgrnas.
    It is used for multiplexing while the n_sgrnas is the amount of guides in a single plasmid (the multiplex) and
    the n_multiplx is the number of groups of sg with the same 'best candidate' to return
    The algorithm takes the best sgRNA and match it with different guides to create m_group of multiplex each one with n_sgrnas guides.
    It returns a BestSgGroup object with the attribute 'subgroups' that store a list of subgroups the length of m_groups,
     each subgroup is a SubgroupRes object with n_sgrnas as the amount of candidates in its candidates_list
    Args:
        subgroup: SubgrouRes object
        n_multiplx: number of groups of guides (the number of BestSgGroups to create)
        n_sgrnas: number of guides in each group (for multiplexing)

    Returns:
            BestSgGroup object
    """
    # make a copy of subgroup
    subgroup_temp = copy.deepcopy(subgroup)
    # randomly shuffle the gRNA list
    random.shuffle(subgroup_temp.candidates_list)
    # if no more candidates left stop the search
    if not subgroup_temp.candidates_list:
        return None
    # get the first group of sgRNAs, the best guide in the group and the positions of the candidate
    first_multiplex = choose_multiplx(subgroup_temp, n_sgrnas)
    # check if a multiplex is found
    if not first_multiplex:
        return None
    # save the result to a BestSgGroup object, this object is design to hold all multiplex of the same 'best' sgRNA
    current_best = BestSgGroup()
    # store the subgroupres object with the candidates
    current_best.subgroups_list = [first_multiplex[0]]
    # add the 'best' and the postion
    current_best.best_candidate, pos_dict = first_multiplex[1], first_multiplex[2]
    # store a list of all candidates
    current_best.all_candidates = copy.copy(current_best.subgroups_list[0].candidates_list)
    while n_multiplx > 1:
        # remove the found sg from the subgroup (recreate it without them)
        subgroup_temp.candidates_list = [can for can in subgroup_temp.candidates_list if
                                         can not in current_best.all_candidates]
        # check that the are candidates left in the subgroup
        if not subgroup_temp.candidates_list:
            return current_best
        # choose the next group of sgRNA that will be joined with the 'best guide' found above
        try:
            # get the SubgroupRes, best_candidate and the updated positions list
            multiplx_group, best_candidate, pos_dict = choose_multiplx(subgroup_temp, n_sgrnas,
                                                                        current_best.best_candidate, pos_dict)
            # add the SubgrouRes object containing the list of guides to the results
            current_best.subgroups_list.append(multiplx_group)
            current_best.all_candidates += [can for can in multiplx_group.candidates_list if
                                            can not in current_best.all_candidates]
            n_multiplx -= 1
        except TypeError:
            # print(f"No more candidate in group {current_best.best_candidate.seq} in node {subgroup_temp.name}")
            return current_best
    return current_best


def add_best_candidate_to_multiplx(output_dict: Dict):
    """
    This function add the candidate object of the 'Best' candidate for each subgroup object
    it is done in order to coolect all multiplx of the same 'Best gRNA' together in later steps.
    Args:
        output_dict: the dictionary of {node_name:{best_seq:BestSgGroup}}

    Returns:

    """
    for node in output_dict.keys():
        for bestgroup in output_dict[node].values():
            for multiplx in bestgroup.subgroups_list:
                for cand in multiplx.candidates_list:
                    cand.best_candidate = bestgroup.best_candidate


def filter_duplicates_from_final_dict(res_dict: Dict) -> Dict:
    """
    This function remove duplicate multiplx (where the sequence of all gRNA are the same)
    It will leave the the multiplx that target the most number of genes (highest internal node on the gene tree)
    And in cases where both target the same genes it will leave the one with the smallest amount of genes in the node
    (the lowest internal node)
    Args:
        output_dict: Dictionary of internal_node: dict of best_seq:BestSgGroup

    Returns:
         A dictionary of (tuple of seqs):SubgroupRes (multiplex)

    """
    multiplx_dict = {}
    # go over the dictionary on nodes and collect the multiplx by there sequence to a new output dict
    for node in res_dict.keys():
        for bestgroup in res_dict[node].values():
            for multiplx in bestgroup.subgroups_list:
                # get the sequences of the gRNAs in the multiplx
                multi_seqs = tuple(sorted([can.seq for can in multiplx.candidates_list]))
                # get the number of genes targeted by the gRNAs
                n_genes_targeted_current = sum([len(can.genes_score_dict.keys()) for can in multiplx.candidates_list])
                # if a multiplx with the same sequnce is in the output dict compare the number of genes each is targeting
                if multi_seqs in multiplx_dict.keys():
                    n_genes_targeted = sum([len(can.genes_score_dict.keys()) for can in multiplx_dict[multi_seqs].candidates_list])
                    # if the current multiplx target more genes replace the existing one in the dictionary
                    if n_genes_targeted_current > n_genes_targeted:
                        multiplx_dict[multi_seqs] = multiplx
                    # if the number of genes targeted is equal, choose the one with the lower intrnal node
                    elif len(multiplx.genes_in_node) < len(multiplx_dict[multi_seqs].genes_in_node):
                        multiplx_dict[multi_seqs] = multiplx
                else:
                    multiplx_dict[multi_seqs] = multiplx
    return multiplx_dict


def select_top_sg_in_node(multiplx_dict: Dict, sg_per_node=4) -> Dict:
    """
    This function take the dictionary of chips results {(seqs_tuple):subgroupres} and output dictionary of
    {internal_node:n top gRNAs}
    First it calculate the total score of each multiplex (sun of cut expectation) and add it as an attribute to each multiplex,
    Than it re-create the tree structure by assigning to each internal node its candidates and finaly it will return
    the top 'sg_per_gene' for each internal node.
    Args:
        multiplx_dict: dictionary of multiplxes
        sg_per_node: the number of top gRNAs to select from each internal node

    Returns: final dictionary results with n best gRNAs per internal node

    """
    res_dict = {}
    # create a list of all multiplexes
    multiplx_lst = list(multiplx_dict.values())
    # shuffle list
    random.shuffle(multiplx_lst)
    # add score of all gRNAs in muktpilx
    for multiplx in multiplx_lst:
        total_score = 0
        for can in multiplx.candidates_list:
            total_score += can.cut_expectation
        multiplx.total_score = total_score

    # rearrange the multiplx to group genes in node
    for multiplx in multiplx_lst:
        if str(multiplx.genes_in_node) not in res_dict:
            res_dict[str(multiplx.genes_in_node)] = [multiplx]
        else:
            res_dict[str(multiplx.genes_in_node)].append(multiplx)

    # take n multiplx from each gene group
    #sort multiplexes by total score
    res_dict = {key:sorted(value, key=lambda x: x.total_score, reverse=True) for key,value in res_dict.items()}
    # output results with the top multiplex as defined in the 'sg_per_node' variable
    res_dict = {key:value[0:sg_per_node] for key, value in res_dict.items()}

    return res_dict


def chips_main(crispys_output_path: str = None,
               crispys_output_name: str = "crispys_output",
               chips_output_name: str = "chips",
               genome_by_chr_path: str = None,
               pam_file_path: str = None,
               gff_file_path: str = None,
               max_number_of_mismatches: int = 4,
               lower_intersect_limit: int = 10,
               upper_intersect_limit: int = 20,
               number_of_groups: int = 20,
               n_multiplx: int = 5,
               n_sgrnas: int = 2,
               threads: int = 1,
               number_of_singletons: int = 5,
               scoring_function: str = "moff",
               restriction_site: str = "None",
               sg_per_node: int = 5,
               include_all_family_gene: bool = False):
    """
    This is the main function of crispys-chips, it read the result of crispys and apply on them the chips algorithm.
    First, an off-target search is done and only gRNAs that pass it will be processed in further steps.
    Than algorithm looks for n (usualy n=2) size group of gRNA that aim to capture all genes in the family input.
    This group of n gRNAs is a multiplx, the process of selecting multiplx is, in short, to first find the gRNA that
    capture the most genes and than to add to it gRNAs the capture genes not targeted by the 'Best' gRNA.
    For example, when n=2 and the family has 4 genes we first select gRNA that capture the most genes in the family,
    lets say that our 'Best' gRNA captured g1, g2 and g3 than we will look for another gRNA (secondary) to capture g4.
    The process of adding different secondary gRNA that will complement the 'Best' gRNA is repeated for 'n_multiplx'
    times (default=5) and than another 'Best' candidate is selected and the whole process is repeated for 'number_of_groups'
    times (default=20)
    Args:
        crispys_output_path: path to crispys output folder
        crispys_output_name: name given to crispys output (set in crispys run)
        chips_output_name: name to give the output
        restriction_site: a sequence of dna recognized by the enzyme used in the transformation reaction, if a gRNA
                        have this sequence it will be exclude from output
        genome_by_chr_path: path to the folder the have the fasta file of reference genome each file holds one chromosome
        pam_file_path: path to folder that have a pan file, this file is used to find off-targets, will be created if not exist.
        gff_file_path: path to gff annotation file, used to find the genomic element of off-targets.
        max_number_of_mismatches: maximum number of mismatches to look for in off-target search.
        lower_intersect_limit: Lower limit for the intersection between the off-target and a genomic region.
        upper_intersect_limit: Upper limit for the intersection between the off-target and a genomic region.
        number_of_groups: The number of groups of 'best' multiplxes to output
        n_multiplx: The number of multiplxes for each 'Best' gRNA
        n_sgrnas: the number of gRNA in multiplx
        threads: number if cpu to use with crispritz
        number_of_singletons: how many singletons to add to each group of gRNAs of internal node
        scoring_function: the scoring function to use in off-target search
        sg_per_node: the number of gRNA to select from each subgroup of genes (internal node)

    Returns:
        write several output file including csv with final multiplex
    """

    # read CRISPys results
    pickle_name = [file for file in os.listdir(crispys_output_path) if file.endswith(f"{crispys_output_name}.p")][0]
    pickle_file = os.path.join(crispys_output_path, pickle_name)
    list_of_subgroup_results = pickle.load(open(pickle_file, 'rb'))
    # Check the results if there is only one gene
    check_if_one_gene(list_of_subgroup_results)
    # Convert the list of internal node results into a list of CandidateWithOffTargets object.
    list_of_candidates = concat_candidates_from_subgroups(list_of_subgroup_results)
    # remove restriction site
    list_of_candidates = remove_candidates_with_restriction_site(list_of_candidates, restriction_site=restriction_site)
    # search off targets
    off_target_search_crispritz_main(list_of_candidates, globals.crispritz_script_path, crispys_output_path,
                                      gff_file_path, genome_by_chr_path,
                                      pam_file_path, max_number_of_mismatches, lower_intersect_limit,
                                      upper_intersect_limit, threads, scoring_function)

    # optionally, get the original family genes name (before the split to smaller families) in order to ignore off-taeget
    # hitting these genes
    family_name = list_of_subgroup_results[0].family_name
    if include_all_family_gene:
        familiy_genes_dict = pickle.load(open(globals.familiy_genes_dict_path, 'rb'))
        genes_to_ignore_set = set(familiy_genes_dict[family_name])
    else:
        genes_to_ignore_set = {}

    # filter by genomic region
    list_of_candidates_filtered = filter_main(list_of_candidates, crispys_output_path, globals.feature_score_dict,
                                               globals.ignore_regions, f"{crispys_output_name}_{chips_output_name}",
                                               genes_to_ignore_set)

    # re-create the list of subgroup (crispys output from the sgRNAs that pass the off-target filtering)
    list_of_subgroup_no_offtargets = recreate_subgroup_lst(list_of_candidates_filtered, list_of_subgroup_results)
    # insert singleton subgroup to subgroups without singleton and create a list of subgroups without sinlgetons
    new_subgroups_lst = add_singletons_to_subgroup(list_of_subgroup_no_offtargets, number_of_singletons)

    # Stop if the list is empty
    if not new_subgroups_lst:
        print("No input for Chips\n")
        quit()
##### start Chips ####################################################################################################
    global pos_dict
    global guide_seq_dict
    output_dict = {}
    # sort the subgroup list according the number of genes in node (from low to high)
    new_subgroups_lst.sort(key=lambda x: len(x.genes_in_node), reverse=True)
    # go over each group of results (for each internal node in gene tree)
    for subgroup in new_subgroups_lst:
        subgroup_dict = subgroup2dict(subgroup)
        # check if no candidate in node result
        if len(subgroup_dict) == 0:
            print(f"No CRISPys results for node {subgroup.name}")
            continue
        # initiate results dictionary
        bestsgroup_dict = {}
        # choose multiplex groups
        n = number_of_groups
        while n > 0:
            # get results for a group of guides with the same 'best' guide
            bestsgroup = get_best_groups(subgroup, n_multiplx, n_sgrnas)
            # if there are no more candidate it will return None
            if not bestsgroup:
                break

            # check if 'best' we have got has the same position as other in the group and if so remove it from the
            # results dictionary and the group list
            # if it is the first best group, create the positional dictionary, otherwise compare with previous positions
            if not bestsgroup_dict:
                pos_dict = check_overlap_positions(bestsgroup.best_candidate)
                # store the group of guides in a dictionary with best_sg_seq: BestSgGroup object
                bestsgroup_dict[bestsgroup.best_candidate.seq] = bestsgroup
            else:
                skip_candidate = check_overlap_positions(bestsgroup.best_candidate, pos_dict)
                if skip_candidate:
                    subgroup.candidates_list.remove(subgroup_dict[bestsgroup.best_candidate.seq])
                    continue
                else:
                    bestsgroup_dict[bestsgroup.best_candidate.seq] = bestsgroup

            # remove the 'best candidate' from the list of candidates
            subgroup.candidates_list.remove(subgroup_dict[bestsgroup.best_candidate.seq])
            if not subgroup.candidates_list:
                break
            n -= 1
        if bestsgroup_dict:
            output_dict[subgroup.name] = bestsgroup_dict

    #add the candiate object for each multiplex subgroupres object
    add_best_candidate_to_multiplx(output_dict)
    # filter duplicate multiplexs
    multiplx_dict = filter_duplicates_from_final_dict(output_dict)
    # get output by selecting minimum of n gRNA for each gene
    final_dict = select_top_sg_in_node(multiplx_dict, sg_per_node)

    # save results to pickle
    with open(f"{os.path.join(crispys_output_path, chips_output_name)}.p", 'wb') as f:
        pickle.dump(final_dict, f)
    # write results to csv file
    create_output_multiplex(crispys_output_path, final_dict, number_of_groups,
                             n_multiplx, n_sgrnas, chips_output_name)



def parse_arguments(parser_obj: argparse.ArgumentParser):
    """
        using a pars_obj object this function parses command line strings into python objects. the chosen parameters for the
        algorithm run are taken from this function.
        :param parser_obj: object for parsing command line strings into python objects
        :return: parsed parameters for the algorithm run
        :rtype: argparse.Namespace
        """

    parser_obj.add_argument('--crispys_output_path', '-crispys', type=str,
    help='The path to CRISPys output folder')

    parser_obj.add_argument('--crispys_output_name', '-crispys_name', type=str,
                             help='crispys output name')

    parser_obj.add_argument('--chips_output_name', '-chips_name', type=str,
                             help='Chips output name')

    parser_obj.add_argument('--genome_by_chr_path', '-genome_chr', type=str,
                             help='The path to the folder with chromosome fasta files (for crispritz)',
                             default=None)

    parser_obj.add_argument('--pam_file_path', '-pam_file', type=str, help='The path to the folder of pam file'
                                                                            ' (for crispritz)', default=None)

    parser_obj.add_argument('--gff_file_path', '-gff', type=str, help='The path to the gff file', default=None)

    parser_obj.add_argument('--max_number_of_mismatches', '-n_mm', type=int, default=4,
                             help='Number of mismatches to allow when searching for off-targets')

    parser_obj.add_argument('--lower_intersect_limit', '-lower', type=int, default=10,
                             help='The minimum overlap between guide and off-target')

    parser_obj.add_argument('--upper_intersect_limit', '-upper', type=int, default=20,
                             help='The maximum overlap between guide and off-target')

    parser_obj.add_argument('--number_of_groups', '-groups', type=int, default=20,
                             help='Number of Best Groups i.e. number of BestSgGroup objects')

    parser_obj.add_argument('--n_multiplx', '-n_multi', type=int, default=5,
                             help="Number of multiplexes for each group of 'Best guide'")

    parser_obj.add_argument('--n_sgrnas', '-sgrnas', type=int, default=2,
                             help='Number of sgRNAs in each multiplex')

    parser_obj.add_argument('--threads', '-th', type=int, default=1,
                             help='Number of threads to use in parallel when using CRISPRitz')

    parser_obj.add_argument('--number_of_singletons', '-n_singletons', type=int, default=5,
                             help='Number of singletons to select from results. '
                                  'for each gene in internal node this number of the best singletons will be added')

    parser_obj.add_argument('--scoring_function', '-scoring', type=str, help="The scoring function to use",
                             default="gold-off", choices=['gold_off', 'moff'])

    parser_obj.add_argument('--restriction_site', '-restriction', type=str,
                         help='The sequence of restriction site')

    parser_obj.add_argument('--sg_per_node', '-node_sg', type=int, help="The minimum number of gRNAs for any gene",
                         default=5)

    arguments = parser_obj.parse_args()
    return arguments


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parse_arguments(parser)
    chips_main( crispys_output_path=args.crispys_output_path,
                crispys_output_name=args.crispys_output_name,
                chips_output_name=args.chips_output_name,
                restriction_site=args.restriction_site,
                genome_by_chr_path=args.genome_by_chr_path,
                pam_file_path=args.pam_file_path,
                gff_file_path=args.gff_file_path,
                max_number_of_mismatches=args.max_number_of_mismatches,
                lower_intersect_limit=args.lower_intersect_limit,
                upper_intersect_limit=args.upper_intersect_limit,
                number_of_groups=args.number_of_groups,
                n_multiplx=args.n_multiplx,
                n_sgrnas=args.n_sgrnas,
                threads=args.threads,
                number_of_singletons=args.number_of_singletons,
                scoring_function=args.scoring_function,
                sg_per_node=args.sg_per_node )