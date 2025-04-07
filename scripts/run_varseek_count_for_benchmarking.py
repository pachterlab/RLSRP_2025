import varseek as vk

_, keyword_args = vk.utils.make_positional_arguments_list_and_keyword_arguments_dict()
fastqs = keyword_args.pop("fastqs", None)
vk.count(fastqs, **keyword_args)