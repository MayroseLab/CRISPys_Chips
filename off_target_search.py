import os
from typing import List, Dict
import subprocess
import time
from os.path import dirname, abspath
import numpy as np
from numpy import clip
import re

from pybedtools import BedTool
from Candidate import Candidate
from crunch_classes import CandidateWithOffTargets, OffTarget

import SubgroupRes
import gold_off
import globals
import random
import pandas as pd


def create_crispritz_input_file(list_of_candidates: List, crispritz_path: str) -> str:
    """
    :param list_of_candidates: A list of CandidateWithOffTargets objects
    :param crispritz_path: A path to the crispys result folder where a folder for crispritz will be created
    :return: A path to the input for cas-offinder and will write the input file to cas_offinder_path
    """
    crispritz_infile = os.path.join(crispritz_path, 'crispritz_infile.txt')
    out = ""
    for candidate in list_of_candidates:  # go over each candidate and get the guide sequence
        out += f"{candidate.seq}NNN\n"
    with open(crispritz_infile, 'w') as f:
        f.write(out)
    return crispritz_infile


def create_pam_file(pam_file_path: str) -> str:
    """
    create the pam file for crispritz (for NGG pam)
    Args:
        pam_file_path: path to folder location

    Returns:
        pam file path
    """
    pam_file = f"{pam_file_path}/pamNGG.txt"
    if not os.path.exists(pam_file):
        # os.makedirs(pam_file_path)
        with open(pam_file, "w") as f:
            f.write("NNNNNNNNNNNNNNNNNNNNNGG 3")
        return pam_file
    else:
        return pam_file


def run_crispritz(list_of_candidates: List, crispritz_script_path: str, crispys_output_path: str,
                  genome_by_chr_path: str,
                  pam_file_path, max_number_of_mismatches: int = 4, threads: int = 1) -> pd.DataFrame:
    """
    This function runs crispritz and returns its output
    nucleotide codes (e.g. 'NRG' matches both NGG and NAG).
    :param crispritz_script_path: path to crispritz python script path.
    :param list_of_candidates: A list of CandidateWithOffTargets objects
    :param crispys_output_path: A path containing the output of CRISPys
    :param genome_by_chr_path: path to the folder where the files of each chromosome fasta file.
    :param pam_file_path: a path to folder where pam file (created if not exist)
    :param max_number_of_mismatches: The maximum number of mismatches allowed between a sgRNA
     and a potential off-target.
    :param threads: number of threads to use
    :return: The output of crispritz as pd.DataFrame, where each row is a potential offtarget.
    """
    crispritz_path = f"{crispys_output_path}/crispritz"
    os.makedirs(crispritz_path, exist_ok=True)
    # create crispritz input file
    crispritz_infile = create_crispritz_input_file(list_of_candidates, crispritz_path)
    pam_file = create_pam_file(pam_file_path)
    # run crispritz
    os.system(
        f"python {crispritz_script_path} search {genome_by_chr_path}/ {pam_file} {crispritz_infile} {crispritz_path}/crispritz -mm {max_number_of_mismatches} -r -th {threads}")
    # get results
    crispritz_results = pd.read_csv(f"{crispritz_path}/crispritz.targets.txt", sep="\t")
    print(f"Number of off-targets: {crispritz_results.shape[0]}")
    return crispritz_results


def get_off_target(x, sequence_to_candidate_dict):
    """
    A function to use with apply on crispritz reault table
    it takes a row of crispritz results and a dictionary of sequence:candidate, and make an OffTarget
    object from the crispritz results and add it to the candiddate offtargets list
    Args:
        x: a row in crispritz resuls table
        sequence_to_candidate_dict: a dictionary of sequence:candidate
    Returns:
        none
    """
    candidate = sequence_to_candidate_dict[x['crRNA'][:20]]
    off_target_seq = x['DNA'].upper()
    matches = {"A", "C", "G", "T"}
    # Check if the off target contains ACGT letters only before creating the OffTarget object
    if all(nucl in matches for nucl in off_target_seq):
        off_target = OffTarget(off_target_seq, x['Chromosome'], int(x['Position']), x['Direction'],
                               int(x['Mismatches']))

        candidate.off_targets_list.append(off_target)
    return


def add_crispritz_off_targets(crispritz_results, sequence_to_candidate_dict: Dict):
    """
    This function adds all found off-targets to each CandidateWithOffTargets using the crispritz results
    :param crispritz_results: The output of crispritz as a pd datatable, where each row is a potential offtarget.
    :param sequence_to_candidate_dict: sequence -> a CandidateWithOffTargets object with the proper sequence
    :return: None
    """
    # apply the 'get_off_target' function on each row in the crispritz table results
    crispritz_results.apply(get_off_target, args=(sequence_to_candidate_dict,), axis=1)
    return


def run_gold_off(sg_list: list, off_targets_list: list):
    """
    This function get 2 lists of sequences, first a list of crispr genomic target (with pam)
    and second list with offtargets sequences we want to calculate relative to the target.
    both lists needs to be with the same length
    Args:
        sg_list: list of sgrna target
        off_targets_list: list of off-targets

    Returns:
        a list of gold-off scores
    """

    xgb_model_path = globals.goldoff_model_path
    list_of_scores = gold_off.predict(sg_list, off_targets_list, xgb_model_path, include_distance_feature=True,
                                      n_process=1, model_type="regression")
    return clip(list_of_scores, 0, 1)


def add_scores_to_off_targets(candidates_list, scoring='gold-off'):
    """
    This function take a list of candidate objects and add the results of the scoring function to each OffTarget object
    inside CandidateWithOffTargets.
    Args:
        candidates_list: list of candidates after offtarget sequence was added
        scoring: name of scoring function algorithm

    Returns:
        add a score for each OffTarget
    """
    if scoring == 'gold-off':
        for candidate in candidates_list:
            off_seq_lst = [off.seq for off in candidate.off_targets_list]
            if not off_seq_lst:
                continue
            sgRNA_lst = [f"{candidate.seq}NGG"] * len(candidate.off_targets_list)
            # get gold off results
            gold_off_scores = run_gold_off(sgRNA_lst, off_seq_lst)
            # add results to OffTarget object
            for off, score in zip(candidate.off_targets_list, gold_off_scores):
                off.score = score

    elif scoring == "moff":
        from MOFF.MoffLoad import load_moff
        from MOFF.MOFF_prediction import MOFF_score
        load_moff()
        for candidate in candidates_list:
            off_seq_lst = [off.seq for off in candidate.off_targets_list]
            if not off_seq_lst:
                continue
            sgRNA_lst = [f"{candidate.seq}NGG"] * len(candidate.off_targets_list)
            # calculate moff scores
            moff_scores = MOFF_score(globals.moff_mtx1, globals.moff_mtx2, sgRNA_lst, off_seq_lst)
            for off, score in zip(candidate.off_targets_list, moff_scores):
                # get the complement score.
                off.score = 1 - score
    elif scoring == "random":
        for candidate in candidates_list:
            off_seq_lst = [off.seq for off in candidate.off_targets_list]
            # get random results
            random_scores = [random.random() for rand in range(len(off_seq_lst))]
            # add results to OffTarget object
            for off, score in zip(candidate.off_targets_list, random_scores):
                off.score = score
    else:
        print("No off-target scoring function selected!")
        exit()

def run_bedtools(list_of_candidates: List, gff_file_path: str,
                 lower_intersect_limit: int = 10, upper_intersect_limit: int = 20) -> BedTool:
    """
    This function assumes that the genome used for the CRISPys input and the annotation file are from
    the same database.
    :param list_of_candidates: A list of CandidateWithOffTargets objects
    :param gff_file_path: A path to a gff file containing structural annotation of a given genome.
    IMPORTANT: make sure that the chromosome naming used in the annotation file and the genome file are the same.
    :param lower_intersect_limit: Lower limit for the intersection between the off-target and a genomic region.
    :param upper_intersect_limit: Upper limit for the intersection between the off-target and a genomic region.
    :return: A BedTool object containing the intersection between the off-targets and the genomic regions
    in the gff file
    """
    off_targets_bed = create_off_targets_bed(list_of_candidates, lower_intersect_limit, upper_intersect_limit)
    annotation_bed = BedTool(gff_file_path)
    intersect_bed = off_targets_bed.intersect(annotation_bed, wb=True, stream=True).saveas()
    return intersect_bed


def add_genomic_regions_to_off_targets(intersect_bed, sequence_to_candidate_dict):
    """
    :param intersect_bed: A BedTool object containing the intersection between the off-targets and the genomic regions
    in the gff file
    :param sequence_to_candidate_dict: sequence -> a CandidateWithOffTargets object with the proper sequence
    :return:
    """
    for interval in intersect_bed:
        sgRNA_candidate = sequence_to_candidate_dict[interval.name]

        # get feature metadate (last column of gff)
        attributes = interval.fields[13]

        # get the off target by the index created in the bed object
        off_target = sgRNA_candidate.off_targets_list[int(interval.score)]
        genomic_region = interval.fields[7]
        # dont write 'chromosome' as genomic region
        if genomic_region.lower() == 'chr' or genomic_region.lower() == 'chromosome':
            continue
        # add the genomic region to offtarget.genomic_region set
        off_target.genomic_regions.add(genomic_region)
        # if the genomic region is a gene, add its name
        if genomic_region.lower() == "gene":
            gene_name = re.search("ID=([^;]*)", attributes, re.IGNORECASE)
            # the gene name will be added only if it can be found
            try:
                off_target.genes_covered.add(gene_name[1].upper())
            except:
                print(f"\ncould not get gene name from gff!\nchange the regular expression in function "
                      f"'add_genomic_regions_to_off_targets' in 'off_target_search.py'\n"
                      f"the attribute of the gff is:\n{attributes}\n")
                break
    return


def create_off_targets_bed(list_of_candidates, lower_intersect_limit: int = 10,
                           upper_intersect_limit: int = 20) -> BedTool:
    """
    :param list_of_candidates: A list of CandidateWithOffTargets objects.
    :param lower_intersect_limit: Lower limit for the intersection between the off-target and a genomic region.
    :param upper_intersect_limit: Upper limit for the intersection between the off-target and a genomic region.
    :return: A BedTool object containing all off-target sequences
    """
    off_targets_bed_string = ""
    for candidate in list_of_candidates:
        for i, off_target in enumerate(candidate.off_targets_list):
            off_targets_bed_string += f"{off_target.chromosome} {off_target.start_position + lower_intersect_limit} {off_target.start_position + upper_intersect_limit} {candidate.seq} {i}\n"
    return BedTool(off_targets_bed_string, from_string=True)


def create_sequence_to_candidate_dict(list_of_candidates):
    """
    This function takes a list of sgRNA candidates, and returns a sequence to candidate dictionary.
    :param list_of_candidates: a list of sgRNA candidates
    :return: a dictionary: sequence -> a CandidateWithOffTargets object where candidate.seq = sequence.
    """
    assert all([isinstance(element, Candidate) for element in list_of_candidates])
    sequence_to_candidate_dict = {}
    for i, candidate in enumerate(list_of_candidates):
        sequence_to_candidate_dict[candidate.seq] = list_of_candidates[i]
    return sequence_to_candidate_dict


def off_target_search_crispritz_main(list_of_candidates: List, crispritz_script_path: str, crispys_output_path: str,
                                     gff_file_path: str, genome_by_chr_path: str, pam_file_path: str,
                                     max_number_of_mismatches: int = 4, lower_intersect_limit: int = 10,
                                     upper_intersect_limit: int = 20, threads: int = 1,
                                     scoring_function: str = "gold-off"):
    """
    CRISPR Off-target search algorithm for sgRNAs given by CRISPys.
    :param genome_by_chr_path: path to the folder where the files of each chromosome fasta file.
    :param list_of_candidates: A list of Candidate objects.
    :param crispys_output_path: A path containing the output of CRISPys.
    of a given gene family.
    :param genome_path: A path to the genome file where off-targets are searched.
    :param gff_file_path: A path to a gff file containing structural annotation of a given genome.
    IMPORTANT: make sure that the chromosome naming used in the annotation file and the genome file are the same.
    :param pam: A string containing the desired protospacer adjacent motif (PAM), using the IUPAC
    nucleotide codes (e.g. 'NRG' matches both NGG and NAG).
    :param max_number_of_mismatches: The maximum number of mismatches allowed between a sgRNA
     and a potential off-target.
    :param lower_intersect_limit: Lower limit for the intersection between the off-target and a genomic region.
    :param upper_intersect_limit: Upper limit for the intersection between the off-target and a genomic region.
    :param threads: number of threads to use in parallel
    :return: Updates the candidates in subgroup_results.
    """
    sequence_to_candidate_dict = create_sequence_to_candidate_dict(list_of_candidates)
    # run crispritz
    crispritz_results = run_crispritz(list_of_candidates, crispritz_script_path, crispys_output_path,
                                      genome_by_chr_path, pam_file_path, max_number_of_mismatches, threads)
    # add off targets to candidate
    add_crispritz_off_targets(crispritz_results, sequence_to_candidate_dict)

    # get off-targets results
    add_scores_to_off_targets(list_of_candidates, scoring=scoring_function)
    intersect_bed = run_bedtools(list_of_candidates, gff_file_path, lower_intersect_limit=lower_intersect_limit,
                                 upper_intersect_limit=upper_intersect_limit)

    add_genomic_regions_to_off_targets(intersect_bed, sequence_to_candidate_dict)


