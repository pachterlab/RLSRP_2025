import varseek as vk

_, keyword_args = vk.utils.make_positional_arguments_list_and_keyword_arguments_dict()
fastqs = keyword_args.pop("fastqs", None)
vk.count(fastqs, **keyword_args)

# vk.count("tmp/reads_1000_fastq.fastq", index="/home/jmrich/Desktop/RLSRWP_2025/data/vk_ref_out/vcrs_index.idx", t2g="/home/jmrich/Desktop/RLSRWP_2025/data/vk_ref_out/vcrs_t2g_filtered.txt", technology="bulk", threads=16, k=51, out="tmp/vk_count_threads_16_reads_1000_out", kb_count_reference_genome_out_dir="tmp/kb_count_reference_genome_out_dir_1000", disable_summarize=True)