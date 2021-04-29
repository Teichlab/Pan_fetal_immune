### Splitting rules ###
splitting_scheme = {
    "MYELOID_LYMPHOID":["MYELOID","LYMPHOID"],
    'LYMPHOID':["NKT", "B"],
    "B":["B_CLEAN"],
    "NKT":["NKT_CLEAN"]
}

progenitor_labels = ['MPP', 'MEP', 'HSC', 'MEMP', 'GMP',
 'MYELOID DC PROGENITOR', 'ERYTHROID/MPP', 'MOP',
 'LMPP', 'CMP', 'MPP MYELOID', 'LYMPHOID PROGENITOR']
    
splitting_labels = {
    'LYMPHOID':['PRO B CELL',
 'B CELL',
 'PRE B CELL',
 'NK',
 'ILC PRECURSOR',
 'ELP',
 'SP T CELL',
 'ILC',
 'NK T',
 'PRE PRO B CELL',
 'TREG',
 'NK PROGENITOR',
 'DP T CELL',
 'DN T CELL',
 'LTI/ILC3'] ,
    'NKT':[
 'NK',
 'ILC PRECURSOR',
 'SP T CELL',
 'ILC',
 'NK T',
 'TREG',
 'NK PROGENITOR',
 'DP T CELL',
 'DN T CELL',
 'LTI/ILC3'] ,
    'B':['B CELL', 'PRE B CELL', 'PRE PRO B CELL','ELP'],
    'MYELOID':['KUPFFER CELL',
 'MONOCYTE/MACROPHAGE',
 'DC',
 'MONOCYTE',
 'NEUTROPHIL-MYELOID PROGENITOR',
 'DC PRECURSOR',
 'MACROPHAGE',
 'NEUTROPHIL',
 'MONOCYTE/DC PRECURSOR',
 'NEUTROPHIL PROGENITOR',
 'PROMONOCYTE',
 'OSTEOCLAST',
 'MYELOCYTE',
 'PROMYELOCYTE',
 'LANGERHAN CELLS',
 'MACROPHAGE/DC PRECURSOR',
 'NEUTROPHIL PRECURSOR',
 'CD16 MYELOID']
}

# def assign2split(label, split_name):
#     if label in splitting_labels['lymphoid_labels'] + progenitor_labels:
#         if split_name == "LYMPHOID":
#             split_name
#         if split_name == "B":
            
#             split_label = "LYMPHOID"
#         elif label in splitting_labels['myeloid_labels'] + progenitor_labels:
#             split_label = "MYELOID"
#         else:
#             split_label = np.nan
#     elif split_name == "LYMPHOID":
#         if label in splitting_labels['Bcell_labels'] + progenitor_labels:
#             split_label = "B"
#         elif label in splitting_labels['NKTcell_labels'] + progenitor_labels:
#             split_label = "NKT"
#         else:
#             split_label = np.nan
#     else:
#         split_label = np.nan    
#     return(split_label)
    
        
        
    