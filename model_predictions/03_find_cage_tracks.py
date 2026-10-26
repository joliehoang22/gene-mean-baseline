"""
03_find_cage_tracks.py

Lists the brain and whole-blood CAGE tracks in the Basenji2/Enformer target
set (targets_human.txt). This is how the 7 brain and 4 blood tracks were
picked; the chosen indices are saved in resources/brain_tracks.txt and
resources/blood_tracks.txt.
"""
import os
import pandas as pd

tracks = pd.read_csv(
    os.path.join(os.environ["PROJECT_DIR"], "targets_human.txt"),
    sep="\t"
)

print(f"Total tracks: {len(tracks)}")
print(f"Columns: {tracks.columns.tolist()}")

# find brain CAGE tracks
brain_cage = tracks[
    tracks["description"].str.contains("CAGE", case=False, na=False) &
    tracks["description"].str.contains("brain|cortex|frontal|temporal|occipital|neuron|cerebr", case=False, na=False)
]

print(f"\nFound {len(brain_cage)} brain CAGE tracks:\n")
print(brain_cage[["index", "description"]].to_string())

print(f"\nTrack indices: {brain_cage['index'].tolist()}")

# Brain tracks kept:
# BRAIN_TRACKS = [4708, 5094, 5261, 4680, 4980, 4769, 5112]
# 4708  CAGE:frontal lobe, adult, pool1
# 5094  CAGE:medial frontal gyrus, adult
# 5261  CAGE:medial frontal gyrus - adult
# 4680  CAGE:brain, adult, pool1           (general brain)
# 4980  CAGE:brain, adult                  (general brain)
# 4769  CAGE:Astrocyte - cerebral cortex   (cell type specific)
# 5112  CAGE:Neurons                       (cell type specific)

# Brain tracks dropped:
# 4960, 4969, 4981   fetal (the LBP participants are adults)
# 5051               adrenal cancer cell line (matched the search pattern)
# 5162, 5222         meninges, CSF (wrong tissue)
# 5279, 5280, 5282   newborn

###### BLOOD
cage = tracks[tracks["description"].str.contains("CAGE", case=False, na=False)]

# candidate blood tracks: unsorted whole blood and unstimulated primary blood cells
keep = cage[cage["description"].str.contains(
    r"\bblood\b|\bleukocyte|\bPBMC\b|peripheral blood mononuclear|"
    r"\bmonocyte|\bneutrophil|\beosinophil|\bbasophil|"
    r"\bT cell|\bB cell|\bNK cell|natural killer|\bCD4\+|\bCD8\+|\bCD14\+|\bCD19\+",
    case=False, na=False, regex=True
)]

bad = keep["description"].str.contains(
    r"cell line|leukemia|lymphoma|myeloma|carcinoma|sarcoma|"
    r"cord blood|fetal|newborn|derived|"
    r"treated|stimulat|infection|untreated|mock|"
    r"vessel|endothelial",
    case=False, na=False, regex=True
)

print(keep[~bad][["index", "description"]].to_string())

# This lists 32 tracks: sorted immune cell types plus 4 whole-blood libraries.
# Blood tracks kept (whole blood only, to match bulk whole-blood RNA-seq):
# 4805  CAGE:Whole blood (ribopure), donation1
# 4806  CAGE:Whole blood (ribopure), donation2
# 5067  CAGE:Whole blood (ribopure), donation3
# 4950  CAGE:blood, adult, pool1

# The kept indices are hardcoded in the Enformer predict_personal_*.py scripts and
# read from brain_tracks.txt / blood_tracks.txt by the Basenji2 libs.
