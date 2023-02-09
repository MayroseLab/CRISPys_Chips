import os
import SubgroupRes
from typing import List, Dict


def sub_tree_display(candidates_lst, f):
    header_row = "#sgRNA index,sgRNA,Score,Genes,Genes score,Target site,#mms,Position,Strand,PAM,Family, Off-targets\n"

    f.write(header_row)
    sgRNA_index = 0
    for candidate in candidates_lst:
        sgRNA_index += 1
        num_of_targets = 0
        for targets in candidate.targets_dict.values():
            num_of_targets += len(targets)
        first_gene = 1
        l = list(candidate.targets_dict.items())
        l.sort(key=lambda item: candidate.genes_score_dict[item[0]], reverse=True)

        for gene, targets in l:
            targets.sort(key=lambda target: len(target[1]))
            seen_sites = dict()
            first_target = 1
            for target in targets:
                if first_target == 1 and first_gene == 1:

                    f.write(str(sgRNA_index) + '.,' + candidate.seq + "," + str(candidate.cut_expectation)[:5])
                    f.write("," + gene)
                    score = str(candidate.genes_score_dict[gene])
                    if len(score) > 5:
                        score = score[:5]

                    f.write("," + score)
                    f.write("," + change_mismatch_to_lowercase(target[0], target[1].keys()))
                    f.write("," + str(len(target[1])))
                    # write position (added by Udi 03/11/2022)
                    f.write("," + str(target[3]))
                    # write strand
                    f.write("," + target[4])
                    # write pam
                    f.write("," + target[2])
                    # write family name
                    f.write("," + candidate.subgroup.family_name)
                    # write off_targets
                    top_off_trgt = return_top_off_targets(candidate)
                    f.write(f',"{str(top_off_trgt).strip("[]")}"')
                    f.write("\n")
                    first_target = 0
                    continue
                if first_target != 1:
                    f.write(str(sgRNA_index) + ".,,,,,")

                    f.write(change_mismatch_to_lowercase(target[0], target[1].keys()))
                    f.write("," + str(len(target[1])))
                    # f.write("," + pos)
                    # write position
                    f.write("," + str(target[3]))
                    # write strand
                    f.write("," + target[4])
                    # write pam
                    f.write("," + target[2])
                    # write family name
                    f.write( "," + candidate.subgroup.family_name )
                    # write off_targets
                    top_off_trgt = return_top_off_targets(candidate)
                    f.write(f',"{str(top_off_trgt).strip("[]")}"')
                    f.write("\n")
                if first_target == 1 and first_gene != 1:
                    f.write(str(sgRNA_index) + ".,,,")
                    score = str(candidate.genes_score_dict[gene])
                    if len(score) > 5:
                        score = score[:5]

                    f.write(gene)
                    f.write("," + score)
                    f.write("," + change_mismatch_to_lowercase(target[0], target[1].keys()))
                    f.write("," + str(len(target[1])))
                    # f.write("," + pos)
                    # write position
                    f.write("," + str(target[3]))
                    # write strand
                    f.write("," + target[4])
                    # write pam
                    f.write("," + target[2])
                    # write family name
                    f.write("," + candidate.subgroup.family_name)
                    # write off_targets
                    top_off_trgt = return_top_off_targets(candidate)
                    f.write(f',"{str(top_off_trgt).strip("[]")}"')
                    f.write("\n")

                    first_target = 0
            first_gene = 0

def return_top_off_targets(candidate, n_targets=3):
    """
    This function will take from gRNA off targets list its top n_targets (by score)
    and return them as a list
    Args:
        candidate: object that holds gRNA data
        n_targets: number of off-targets to return

    Returns:

    """
    # create a list of off-targets excluding the ones that target family genes
    off_lst = [off for off in candidate.off_targets_list if not off.in_family_gene]
    # sort the list by off-target score
    sorted(off_lst, key=lambda x: x.score, reverse=True)
    return off_lst[0:n_targets]

def change_mismatch_to_lowercase(target_str, mm_lst):
    """

    :param target_str:
    :param mm_lst:
    :return:
    """
    target_in_lst = list(target_str)
    for place in mm_lst:
        target_in_lst[place] = target_in_lst[place].lower()
    return ''.join(target_in_lst)


def tree_display(path: str, subgroups_lst: list,
                 output_name: str = "crispys_output"):
    """
    This function takes the results of crispys and write the crispys results in a CSV output to the output folder
    Args:
        path: path to output folder
        subgroups_lst: list os subgroup objects
        genes_list: a list of gene sequence
        targets_genes_dict: a dictionary with the all targets and the genes they capture
        output_name: The name for the output of CRISPys. The name of the file would be {crispys_output_name}.csv

    Returns:
    """

    filepath = os.path.join(path, f"{output_name}.csv")
    f = open(filepath, 'w')
    f.write(
        "The designed sgRNAs for the genes in your input are listed in the table below. Every section of the table corresponds to a homologous genes subgroup as specified by the internal nodes of the constructed genes tree.<br>The name of the subgroup and the list of genes are given in the header of each section.\n")
    for subgroup_item in subgroups_lst:
        # create the main table
        sub_tree_display(subgroup_item.candidates_list, f)
    f.close()


def create_output_multiplex_from_bestgroup(path: str, multiplex_dict: Dict, number_of_groups: int,
                            n_with_best_guide: int, n_sgrnas: int, output_name: str):

    """
    This function is used to write the output of multiplex
    Args:
        path: path to output folder
        crispys_res: the 'traditional' crispys output (list of subGroupRes)
        multiplex_dict: the output of crispys-chips. a dictionary of node_name:dictionary of best_sg_seq:BestSgGroup
        output_name:
    Returns:

    """

    filepath = f"{path}/{output_name}.csv"

    f = open(filepath, 'w')
    f.write(f"Off-target filtered results of multiplex run with {number_of_groups} 'Best' groups each one with {n_with_best_guide} "
            f"gRNA and {n_sgrnas} in each multiplex\n")
    # go over each internal node
    for node in multiplex_dict:
        f.write(f"Node:,{node},")
        # get node genes
        genes_set = get_genes_of_nodes(multiplex_dict, node)
        # write the node name
        f.write(f"genes in node:,{str(multiplex_dict[node][list(multiplex_dict[node].keys())[0]].subgroups[0].genes_in_node).strip('[]')}\n")
        f.write(f"genes captured:,{str(genes_set).strip('{}')}\n")
        # go over each 'best guide' group
        for bestseq in multiplex_dict[node].values():
            # write the sequence of best guide
            f.write(f"Group of:,{bestseq.best_candidate.seq}\n")
            # go over each pair (or more) of multiplex and write it to the file
            for subgroup in bestseq.subgroups:
                sub_tree_display(subgroup.candidates_list, f)
        f.write("\n")
    f.close()

def get_genes_of_nodes(multiplex_dict, node):
    """
    This function is used to get the genes that are covered in an chips-crunch output 'node'
    Args:
        multiplex_dict: result of chips-crunch
        node: a key in the chips-crunch results

    Returns:
            set of genes covered by gRNA in the chips-crunch node
    """
    genes = []
    for best in multiplex_dict[node].values():
        for can in best.all_candidates:
            genes.extend(can.genes_score_dict.keys())
    return set(genes)


def write_library_csv(path: str, res_dict: Dict, output_name: str, n_families: int, n_singletons: int,
                      n_no_input: int, families_with_no_output: int):
    """
    This function is used to create the final chips library csv output
    Args:
        path: path to write the file
        res_dict: chips results, a dictionary of best_seq:BestSgGroup object
        output_name: name of output file
        n_families: nu,ber of families that used for CRISPys
        n_singletons: number of families with no chips output because they had only singletons gRNAs
        n_no_input: number of families with no chips output because there was no input gRNAs for chips

    Returns:

    """

    filepath = f"{path}/{output_name}.csv"

    f = open(filepath, 'w')
    f.write(f"# CRISPys-Chips output. out of {n_families} families {n_singletons} had only singletons and {n_no_input}"
            f" had no input for Chips and {families_with_no_output} with no chips/crispys output\n")

    for subgroup_list in res_dict.values():
        # write the sequence of best guide
        # f.write(f"#Group of:,{bestgroup.best_candidate.seq}\n")
        # go over each pair (or more) of multiplex and write it to the file
        for subgroup in subgroup_list:
            sub_tree_display(subgroup.candidates_list, f)
        f.write("\n")
    f.close()


def create_output_multiplex(path: str, chips_res_dict: Dict, number_of_groups: int,
                            n_with_best_guide: int, n_sgrnas: int, output_name: str):

    """
    This function is used to write the output of multiplex
    Args:
        path: path to output folder
        crispys_res: the 'traditional' crispys output (list of subGroupRes)
        multiplex_dict: the output of crispys-chips. a dictionary of node_name:dictionary of best_sg_seq:BestSgGroup
        output_name:
    Returns:

    """

    filepath = f"{path}/{output_name}.csv"

    f = open(filepath, 'w')
    f.write(f"Off-target filtered results of multiplex run with {number_of_groups} 'Best' groups each one with {n_with_best_guide} "
            f"gRNA and {n_sgrnas} in each multiplex\n")
    # go over each internal node
    for node in chips_res_dict:
        f.write(f"Node:,{node}\n")
        # get node genes
        # genes_set = get_genes_of_nodes(multiplex_dict, node)
        # # write the node name
        # f.write(f"genes in node:,{str(multiplex_dict[node][list(multiplex_dict[node].keys())[0]].subgroups[0].genes_in_node).strip('[]')}\n")
        # f.write(f"genes captured:,{str(genes_set).strip('{}')}\n")
        # go over each 'best guide' group
        for subgroup in chips_res_dict[node]:
                sub_tree_display(subgroup.candidates_list, f)
        f.write("\n")
    f.close()


if __name__ == "__main__":
    tree_display("/groups/itay_mayrose/galhyams/1516893877")
