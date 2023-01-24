import copy
from typing import List, Dict


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
    # get a list of all genes targeted
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
    # copy the current gene:n_guides dictionary
    genes_scores = genes_sgcount_dict.copy()
    # count the total number of guides per gene starting with the existing one and adding the current multiplex
    for candidate in multiplx.candidates_list:
        for gene in genes_sgcount_dict.keys():
            if gene in candidate.genes_score_dict.keys():
                genes_scores[gene] += 1
    # check if all gene are above the maximum number of guide, if so return None
    if all([score > max_sg_per_gene for score in genes_scores.values()]):
        return
    # if the maxumum has not been reached, add teh number tio the dictionary and return the multiplx
    else:
        for candidate in multiplx.candidates_list:
            for gene in candidate.genes_score_dict.keys():
                genes_sgcount_dict[gene] += 1
        return multiplx
