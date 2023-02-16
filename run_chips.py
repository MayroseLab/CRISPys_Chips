from chips_main import chips_main

chips_main(crispys_output_path="/groups/itay_mayrose/udiland/crispys_chips_arabidopsis/families/HOM05D000028_11",
            crispys_output_name="moff_0.15",
            chips_output_name="chips_test",
            restriction_site="GGTCTC",
            genome_by_chr_path="/groups/itay_mayrose/udiland/crispys_chips_arabidopsis/ath_split_files",
            pam_file_path="/groups/itay_mayrose/udiland/crispys_chips_arabidopsis/pam",
            gff_file_path="/groups/itay_mayrose/udiland/crispys_chips_arabidopsis/annotation.all_transcripts.all_features_modified.ath.gff3",
            max_number_of_mismatches=4,
            lower_intersect_limit=10,
            upper_intersect_limit=20,
            number_of_groups=20,
            n_multiplx=5,
            n_sgrnas=2,
            threads=1,
            number_of_singletons=5,
            scoring_function="moff",
            sg_per_node=3)


# 7 genes family : /groups/itay_mayrose/udiland/crispys_chips_arabidopsis/families/HOM05D000028_11

# 4 genes familiy: /groups/itay_mayrose/udiland/crispys_chips_arabidopsis/families/HOM05D000028_2

# 3 genes family: /groups/itay_mayrose/udiland/crispys_chips_arabidopsis/families/HOM05D000028_8 "crispys_test"

# 2 genes family: /groups/itay_mayrose/udiland/crispys_chips_arabidopsis/families/HOM05D000177_1