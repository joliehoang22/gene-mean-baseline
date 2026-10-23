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
#       index                                                  description
# 4758   4758                                        CAGE:CD14+ Monocytes,
# 4759   4759                                           CAGE:CD4+ T Cells,
# 4760   4760                                           CAGE:CD8+ T Cells,
# 4764   4764                                   CAGE:Natural Killer Cells,
# 4765   4765                     CAGE:Peripheral Blood Mononuclear Cells,
# 4767   4767                                            CAGE:Neutrophils,
# 4805   4805                     CAGE:Whole blood (ribopure), , donation1
# 4806   4806                     CAGE:Whole blood (ribopure), , donation2
# 4950   4950                                     CAGE:blood, adult, pool1
# 5067   5067                     CAGE:Whole blood (ribopure), , donation3
# 5081   5081                CAGE:CD19+ B Cells (pluriselect), , donation2
# 5082   5082                 CAGE:CD8+ T Cells (pluriselect), , donation1
# 5083   5083                CAGE:CD19+ B Cells (pluriselect), , donation1
# 5084   5084                 CAGE:CD8+ T Cells (pluriselect), , donation2
# 5085   5085                 CAGE:CD8+ T Cells (pluriselect), , donation3
# 5086   5086                CAGE:CD19+ B Cells (pluriselect), , donation3
# 5117   5117                                          CAGE:CD19+ B Cells,
# 5135   5135                                         CAGE:neutrophil PMN,
# 5138   5138                                              CAGE:Basophils,
# 5139   5139                                            CAGE:Eosinophils,
# 5209   5209             CAGE:CD4+CD25+CD45RA- memory regulatory T cells,
# 5210   5210   CAGE:CD4+CD25-CD45RA+ naive conventional T cells expanded,
# 5211   5211     CAGE:CD4+CD25+CD45RA+ naive regulatory T cells expanded,
# 5212   5212    CAGE:CD4+CD25+CD45RA- memory regulatory T cells expanded,
# 5213   5213            CAGE:CD4+CD25-CD45RA+ naive conventional T cells,
# 5214   5214                                   CAGE:CD14-CD16+ Monocytes,
# 5215   5215                                   CAGE:CD14+CD16+ Monocytes,
# 5216   5216  CAGE:CD4+CD25-CD45RA- memory conventional T cells expanded,
# 5217   5217                                   CAGE:CD14+CD16- Monocytes,
# 5218   5218              CAGE:CD4+CD25+CD45RA+ naive regulatory T cells,
# 5219   5219           CAGE:CD4+CD25-CD45RA- memory conventional T cells,
# 5275   5275                           CAGE:gamma delta positive T cells,

# Blood tracks kept (whole blood only):
BLOOD_TRACKS = [4805, 4806, 5067, 4950]
BLOOD_TRACK_NAMES = [
    'CAGE_whole_blood_ribopure_donation1',
    'CAGE_whole_blood_ribopure_donation2',
    'CAGE_whole_blood_ribopure_donation3',
    'CAGE_blood_adult_pool1',
]

# these indices are hardcoded in the Enformer predict_personal_*.py scripts and
# read from brain_tracks.txt / blood_tracks.txt by the Basenji2 libs
