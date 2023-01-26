import pickle
import os
import re
from typing import List, Dict
from chips_main import BestSgGroup
from make_tree_display_CSV import write_library_csv

def get_crispys_families(main_families_path: str, log_file: str):
    """
    This function get the families that CRISpys had been run on
    it reads the log file of crispys and extract a list of families names
    Args:
        main_families_path: the main families folder
        log_file: name of crispys log file, it need to be located in the main families folder

    Returns:
        a list of families names
    """
    with open(os.path.join(main_families_path, log_file)) as f:
        all = f.readlines()
    crispys_families = [re.search(":(.*)", x).group(1).strip() for x in all]
    return crispys_families

def check_if_singletons(family_path: str):
    """
       This function check if chips input had only singletons so it cannot preform the chips algorithm.
       it reads the OU file and look for the words "Only singletons"
       it returns the output of the regular expression search on the OU file
       Args:
           family_path: path to family folder

       Returns: the results of a re.findall it can be an empty list if no match or a list of ["Only singletons"]

       """
    files = os.listdir(family_path)
    # get OU file
    file = str([file for file in files if file.endswith("chips.OU")][0])
    with open(os.path.join(family_path, file)) as f:
        all = f.read()
        singletons = re.findall("Only singletons", all)
        return singletons

def check_chips_input(family_path: str):
    """
    This function check if chips run had no input for chips procedure, that is, there were input for off-target search
    but after it the list of gRNA was empty. it reads the OU file and look for the words "No input"
    it returns the output of the regular expression search on the OU file
    Args:
        family_path: path to family folder

    Returns: the results of a re.findall it can be an empty list if no match or a list of ["No input"]

    """
    files = os.listdir(family_path)
    # get OU file
    file = str([file for file in files if file.endswith("chips.OU")][0])
    with open(os.path.join(family_path, file)) as f:
        all = f.read()
        no_input = re.findall("No input", all)
        return no_input

def check_duplicates(res_dict: Dict):
    subgroup_dict = {}
    for bestgroup in res_dict.values():
        for subgroup in bestgroup.subgroups:
            seqs = tuple([can.seq for can in subgroup.candidates_list])
            if seqs not in subgroup_dict:
                subgroup_dict[seqs] = subgroup
            else:
                print("Du Du Duplicate!")

def write_guides_per_gene(res_dict: Dict, library_output_path: str) -> Dict:
    """
    This function calculate the number of guides targeting each gene and write it to a dictionary of
    gene:number of guides, it than write ti to a pickle file
    Args:
        res_dict: the library reults dictionary
        library_output_path: path to output folder to write the output as pickle file

    Returns:

    """
    gene_nguides_dict = {}
    for bestgroup in res_dict.values():
        for subgroup in bestgroup.subgroups:
            for cand in subgroup.candidates_list:
                genes = cand.genes_score_dict.keys()
                for gene in genes:
                    if gene not in gene_nguides_dict:
                        gene_nguides_dict[gene] = 1
                    else:
                        gene_nguides_dict[gene] += 1
    # write the dictionary to file
    with open(f"{os.path.join(library_output_path, 'guides_per_gene.p')}", "wb") as f:
        pickle.dump(gene_nguides_dict, f)

def write_ngenes_targeted_per_guide(res_dict, library_output_path):
    guide_ngenes_dict = {}
    for bestgroup in res_dict.values():
        for subgroup in bestgroup.subgroups:
            for cand in subgroup.candidates_list:
                genes = tuple(cand.genes_score_dict)
                if cand.seq not in guide_ngenes_dict:
                    guide_ngenes_dict[cand.seq] = genes
                else:
                    guide_ngenes_dict[cand.seq] += genes

    guide_ngens = {guide:len(set(genes)) for guide, genes in guide_ngenes_dict.items()}
    # write the dictionary to file
    with open(f"{os.path.join(library_output_path, 'genes_per_guide.p')}", "wb") as f:
        pickle.dump(guide_ngens, f)

def write_ngenes_targeted_per_multiplx(res_dict, library_output_path):
    multplx_ngenes_dict = {}
    for bestgroup in res_dict.values():
        for subgroup in bestgroup.subgroups:
            n_genes = 0
            for cand in subgroup.candidates_list:
                n_genes += len(cand.genes_score_dict)
            # take the sequences of the multiplx and add the to result dict with the number of genes targeted
            seqs = tuple([can.seq for can in subgroup.candidates_list])
            multplx_ngenes_dict[seqs] = n_genes
    # write the dictionary to file
    with open(f"{os.path.join(library_output_path, 'genes_per_multiplex.p')}", "wb") as f:
        pickle.dump(multplx_ngenes_dict, f)



def collect_fam_res(main_families_path: str, chips_name: str, log_file: str, library_output_path: str):
    """
    The main function that create the Chips gRNA library out of all families outputs
    Args:
        main_families_path: path to main families folder
        chips_name: name of chips output (as set in chips 'run_chips_on_families.py')
        log_file: name of log file created in crispys run and located in the main_families_folder
        library_output_path: where to write the output

    Returns:

    """
    # get from the log file a list of families that crispys was run onto
    crispys_families = get_crispys_families(main_families_path, log_file)
    # get all folders names in main families path
    families = os.listdir(main_families_path)
    families = [fam for fam in families if fam.startswith("HOM")]
    # initiate results dictionary
    final_res_dict = {}
    # initiate counts of missing outputs
    families_with_singletons = 0
    families_with_no_chips_input = 0
    families_with_no_output = 0
    # go over each folder and collect chips results
    for family in families:
        if family not in crispys_families:
            continue
        family_path = os.path.join(main_families_path, family)
        if os.path.isdir(family_path):
            try:
                with open(os.path.join(family_path, f"{chips_name}.p"), 'rb') as f:
                    fam_res_dict = pickle.load(f)

            except FileNotFoundError:
                if check_if_singletons(family_path):
                    families_with_singletons += 1
                    print(f"Only singletons in {family}")
                    continue
                if check_chips_input(family_path):
                    families_with_no_chips_input += 1
                    print(f"No input for chips in {family}")
                    continue
                else:
                    families_with_no_output += 1
                    print(f"No Chips output for {family}\n")
                    continue

        for best_seq in fam_res_dict:
            if best_seq not in final_res_dict:
                final_res_dict[best_seq] = fam_res_dict[best_seq]
            else:
                # print(f"Same 'Best' seq in {family}\n")
                final_res_dict[best_seq].subgroups.extend(fam_res_dict[best_seq].subgroups)
    #print number of multpilx
    n_multiplx = sum([len(best.subgroups) for best in final_res_dict.values()])
    print(f"Total number of multiplexes is {n_multiplx} with {2*n_multiplx} guides")

    check_duplicates(final_res_dict)
    # report the number of guides per gene
    write_guides_per_gene(final_res_dict, library_output_path)
    # report the number of genes targeted by each guide
    write_ngenes_targeted_per_guide(final_res_dict, library_output_path)
    # report the number of genes targeted by each multiplex
    write_ngenes_targeted_per_multiplx(final_res_dict, library_output_path)
    # write to csv
    write_library_csv(library_output_path, final_res_dict, "Chips_library_output", len(crispys_families),
                      families_with_singletons, families_with_no_chips_input, families_with_no_output)

collect_fam_res("/groups/itay_mayrose/udiland/crispys_chips_arabidopsis/families", "chips_moff0.15", "crispys_log.txt","/groups/itay_mayrose/udiland/crispys_chips_arabidopsis")