import os
import globals
import pickle
from utils import createHeaderJob

def check_if_results_exist(fam_path, crispys_output_name, chips_output_name):
    """
    This function checks if a CRISPys output file exists for a particular family.
    Args:
        chips_output_name: The output name for chips
        fam_path: The path to the family directory
        crispys_output_name: The output name for CRISPys
    Returns:

    """
    pickles = [file for file in os.listdir(fam_path) if file.endswith(f"{crispys_output_name}.p")]
    if not pickles:
        return False
    pickle_file = os.path.join(fam_path, pickles[0])
    res = pickle.load(open(pickle_file, "rb"))
    if not res:
        print(f"empty CRISPys results for {pickle_file}")
        with open(os.path.join(fam_path, f"{crispys_output_name}_{chips_output_name}_number_of_sgRNAs.txt"), 'w') as f:
            f.write("0\n" * 6)
        return False
    return True




def run(code_path: str, main_folder_path: str, crispys_name: str = "crispys_output", chips_name: str = "chips",
        n_mm: int = 4, lower: int = 10, upper: int = 20, n_groups: int = 20, n_best: int = 5, n_grna: int = 2,
        n_singletons: int = 5, scoring: str = 'moff', ncpu: int = 4, mem: int = 16, queue: str = "itaym",
        restriction_site: str = "None", minimum_sg: int = 4, chips_output_name: str = "chips"):
    """
    A wrapper function to run CRUNCH on the cluster for multiple folders
    Args:
        chips_output_name:
        crispys_output_name:
        n_mm: number of mismatches for off-target search
        lower: The minimum overlap between guide and off-target
        upper: The maximum overlap between guide and off-target
        overlap: The length of overlap to allow between 2 guides
        best: The final number of sgRNAs that would be picked per internal node
        scoring: The scoring function to use
        ncpu: The number of cores to use in each job
        mem: amount of memory for each job
        queue: Name of queue
        code_path: The path to the CRUNCH code

    Returns:
        Runs CRISPys-Chips on each folder
    """
    families = os.listdir(main_folder_path)
    for family in families:
        fam_path = os.path.join(main_folder_path, family)
        if os.path.isdir(fam_path) and not family.startswith("."):
            if not check_if_results_exist(fam_path, crispys_name, chips_output_name):
                continue
            header = createHeaderJob(fam_path, job_name=f"{family}_{crispys_name}_{chips_output_name}",
                                      ncpu=ncpu, mem=mem,
                                      queue=queue)

            command = f"cd {fam_path}\npython {code_path}/chips_main.py -crispys {fam_path} " \
                      f"-crispys_name {crispys_name} -chips_name {chips_name} -genome_chr {globals.genome_by_chr_path} " \
                      f"-pam_file {globals.pam_file_path} -gff {globals.gff_file} -n_mm {n_mm} -lower {lower}" \
                      f" -upper {upper} -groups {n_groups} -n_guide {n_best} -sgrnas {n_grna} -th {ncpu} -n_singletons " \
                      f"{n_singletons} -scoring {scoring} -restriction {restriction_site} -min_sg {minimum_sg}"
            sh_file_path = os.path.join(fam_path, f"{family}_{crispys_name}_{chips_output_name}.sh")
            with open(sh_file_path, "w") as f:
                f.write(f"{header}\n{command}")
            os.system(f"qsub {sh_file_path}")


if __name__ == '__main__':
    run(code_path="/groups/itay_mayrose/udiland/CRISPys_Chips",
         main_folder_path="/groups/itay_mayrose/udiland/crispys_chips_arabidopsis/families", ncpu=1, queue="itaym",
         crispys_name="moff_0.15", chips_name="chips_moff0.15_1sg", n_mm=4, n_groups=20, n_best=5, n_grna=2, n_singletons=5,
         scoring="moff", restriction_site="GGTCTC", minimum_sg=1)
