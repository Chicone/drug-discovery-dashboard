PRESETS = {
    "m3_popc_build": {
        "environment": "membrane",
        "lipids": {"POPC": 1.0},
        "steps": ["martinize", "insane", "patch_top", "fix_ions"],
    },
    "m3_popc_em": {
        "environment": "membrane",
        "lipids": {"POPC": 1.0},
        "steps": ["martinize", "insane", "patch_top", "fix_ions", "em"],
    },
    "m3_popc_eq": {
        "environment": "membrane",
        "lipids": {"POPC": 1.0},
        "steps": ["martinize", "insane", "patch_top", "fix_ions", "em", "nvt", "npt"],
    },
    "m3_popc_prod_50ns": {
        "environment": "membrane",
        "lipids": {"POPC": 1.0},
        "steps": ["martinize", "insane", "patch_top", "fix_ions", "em", "nvt", "npt",
                  "md"],
        "md_ns": 50,
    },
    "m3_popc_prod_200ns": {
        "environment": "membrane",
        "lipids": {"POPC": 1.0},
        "steps": ["martinize", "insane", "patch_top", "fix_ions", "em", "nvt", "npt",
                  "md"],
        "md_ns": 200,
    },
    "m3_popc_chol_prod_200ns": {
        "environment": "membrane",
        "lipids": {"POPC": 0.7, "CHOL": 0.3},
        "steps": ["martinize", "insane", "patch_top", "fix_ions", "em", "nvt", "npt",
                  "md"],
        "md_ns": 200,
    },
}
