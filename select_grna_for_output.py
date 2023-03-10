import copy
from typing import List, Dict

def recreate_subgroup_lst(list_of_candidates_filtered: List) -> List:
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
    return subgroup_lst

def get_final_multiplx_list(multiplx_dict: Dict, max_sg_per_gene=4) -> List:
    """
    This function filter the results of chips to limit the number of gRNAs per gene such that for the groups of genes in
    the family, the gene with the minimum number of gRNAS will not be above predefined number (max_sg_per_gene)
    Args:
        multiplx_dict: a dictionary of (multplx seqs):SubgroupRes
        max_sg_per_gene: the maximum number of gRNA

    Returns: A list of multiplxes (SubgroupRes objects) after filtering with max_sg_per_gene

    """
    # get list of all multiplx
    multiplx_lst = list(multiplx_dict.values())
    # get a dictionary of all genes and the count of how many times they are tarheted (initiate at 0)
    all_genes = [list(can.genes_score_dict.keys()) for group in multiplx_lst for can in group.candidates_list]
    genes_sgcount_dict = {gene: 0 for genes in all_genes for gene in genes}

##### get all multiplexs as one list sorted by combined score
    # add the total score to each multiplex (sum of cut expectation)
    for multiplx in multiplx_lst:
        total_score = 0
        for can in multiplx.candidates_list:
            total_score += can.cut_expectation
        multiplx.total_score = total_score

    # sort results by total_score
    multiplx_lst.sort(key=lambda x: x.total_score, reverse=True)
### go over the sorted list and collect result such that the minimum number of guide is no more than 'max_sg_per_gene'
    # save the first (highest score) multiplx to results
    output_multiplx_list = [multiplx_lst[0]]
    # update genes count dictionary with the count gene of the first multiplx
    genes_in_selected_sggroup = [gene for cand in output_multiplx_list[0].candidates_list for gene in cand.genes_score_dict.keys()]
    for gene in genes_in_selected_sggroup:
        genes_sgcount_dict[gene] += 1

    # go over the rest of the list and collect sggroups until the minimum nuber of gRNA per genes is equal to max_sg_per_guide
    for multiplx in multiplx_lst[1:]:
        if not update_genes_nsg_dict(multiplx, genes_sgcount_dict, max_sg_per_gene):
            continue
        else:
            output_multiplx_list.append(multiplx)
        # if all genes covered to the desire number of gRNA return the results
        if all([gene_score >= max_sg_per_gene for gene_score in genes_sgcount_dict.values()]):
            return output_multiplx_list
    return output_multiplx_list

def update_genes_nsg_dict(multiplx, genes_sgcount_dict: Dict, max_sg_per_gene: int):
    """
    This function count the number sg for each gene in the multiplex and if not all count are above the threshod for
    number of guide per gene (max_sg_per_gene) it return an updated dictinary of gene: n_sgRNA, on the other hand,
    if all genes have more than the minimum guides per gene it return None
    Args:
        multiplx: a SubgroupRes object holding all sgRNA in a multiplx
        genes_sgcount_dict: a dictionary of each gene and the number of guide targeting it
        max_sg_per_gene: the maximum number of guides allowed for the least targeted gene

    Returns:
        if the minimum has not been achieved, returns the multiplx otherwise return None
    """
    # get the genes that needs to be targeted (they have been targeted in previuos multiplex less than max_sg_per_gene)
    genes2target = [gene for gene in genes_sgcount_dict.keys() if genes_sgcount_dict[gene] < max_sg_per_gene]
    # create a list that will hold the genes that are targeted by current multiplex
    genes_targeted = [gene for candidate in multiplx.candidates_list for gene in candidate.genes_score_dict.keys()]
   # check if
    if not [gene for gene in genes2target if gene in genes_targeted]:
        return
    # if the maximum has not been reached, add the number of hits per gene to the dictionary and return the multiplx
    else:
        for candidate in multiplx.candidates_list:
            for gene in candidate.genes_score_dict.keys():
                genes_sgcount_dict[gene] += 1
        return multiplx

def select_top_sg_in_node(multiplx_dict: Dict, sg_per_node=4) -> List:
    res_dict = {}
    # create alist of all multiplexes
    multiplx_lst = list(multiplx_dict.values())
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
