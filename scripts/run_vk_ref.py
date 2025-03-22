import os
import varseek as vk

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(script_dir), "data")
out = os.path.join(data_dir, "vk_ref_out")
reference_out_dir = os.path.join(data_dir, "reference")
w=47
k=51

vk.ref(
    variants="cosmic_cmc",
    sequences="cdna",
    w=w,
    k=k,
    out=out,
    reference_out_dir=reference_out_dir,
    gtf=True,  # just so that gtf information gets merged into cosmic df
    save_logs=True,
    verbose=True,
    save_variants_updated_csv=True,
    save_variants_updated_exploded_vk_info_csv=True,
    save_variants_updated_filtered_csvs=True,
    dlist_reference_source="t2t",
)