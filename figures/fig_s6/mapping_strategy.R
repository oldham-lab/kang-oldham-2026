library(data.table)

# Sea DFC
file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad")
# Sea DFC metadata
file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/RNAseq/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv")

# Sea MTG
file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
# Sea MTG metadata
file.path(Sys.getenv("DATA_DIR", "/mnt/bdata/gugene"), "datasets/SN_RNAseq/sea_ad_2024/mtg/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13_metadata.csv")

# liu all
"/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025.h5ad"
# liu metadata
"/home/gugene/bdata/@shared/scsn.expr_data/human_expr/postnatal/MIT_AD_Multiomic_Multiregion/snRNA_Matrix.2263395_Cells_July7_2025_sif.csv"

##########
# DFC matching
#############

# Unique SEA DFC subclasses:
#  [1] "Oligodendrocyte" "L6b"             "Pvalb"           "Vip"            
#  [5] "OPC"             "Lamp5 Lhx6"      "Sst"             "L6 CT"          
#  [9] "L6 IT"           "Astrocyte"       "L6 IT Car3"      "L5 IT"          
# [13] "L2/3 IT"         "Microglia-PVM"   "L4 IT"           "Lamp5"          
# [17] "L5/6 NP"         "Pax6"            "Sncg"            "Chandelier"     
# [21] "L5 ET"           "Endothelial"     "Sst Chodl"       "VLMC"   


# Unique Liu PFC subclasses:
#  [1] "Exc L3-5 IT"      "Mic"              "Oli"              "Exc L6 IT"       
#  [5] "Exc L3-4 IT"      "Exc L4-5 IT-2"    "OPC"              "Inh SST"         
#  [9] "Exc L6 CT"        "Inh LAMP5"        "Exc L2-3 IT"      "Inh PVALB"       
# [13] "Ast"              "Inh PAX6"         "Inh VIP"          "Exc L5/6 NP"     
# [17] "Exc L5-6 IT"      "Exc L6b"          "Per"              "Exc L5/6 IT Car3"
# [21] "Exc L5 ET"        "Exc L4-5 IT-1"    "End"              "Fib"             
# [25] "SMC"              "T"                "Exc EC"           "Exc HC"   

# SEA subclasses with Liu counterpart:
#  [1] "Oligodendrocyte" "L6b"             "Pvalb"           "Vip"            
#  [5] "OPC"             "Lamp5 Lhx6"      "Sst"             "L6 CT"          
#  [9] "L6 IT"           "Astrocyte"       "L6 IT Car3"      "L5 IT"          
# [13] "L2/3 IT"         "Microglia-PVM"   "L4 IT"           "Lamp5"          
# [17] "L5/6 NP"         "Pax6"            "Sncg"            "Chandelier"     
# [21] "L5 ET"           "Endothelial"     "Sst Chodl"      

# Unique Liu PFC subclasses:
#  [1] "Exc L3-5 IT"      "Mic"              "Oli"              "Exc L6 IT"       
#  [5] "Exc L3-4 IT"      "Exc L4-5 IT-2"    "OPC"              "Inh SST"         
#  [9] "Exc L6 CT"        "Inh LAMP5"        "Exc L2-3 IT"      "Inh PVALB"       
# [13] "Ast"              "Inh PAX6"         "Inh VIP"          "Exc L5/6 NP"     
# [17] "Exc L5-6 IT"      "Exc L6b"          "Exc L5/6 IT Car3"
# [21] "Exc L5 ET"        "Exc L4-5 IT-1"    "End"                   




##########
# MTG matching
#############

# Unique SEA MTG subclasses:
#  [1] "Oligodendrocyte" "L5 IT"           "L2/3 IT"         "L5/6 NP"        
#  [5] "L4 IT"           "OPC"             "Pvalb"           "Sst"            
#  [9] "Lamp5 Lhx6"      "L6 IT"           "Astrocyte"       "Vip"            
# [13] "VLMC"            "Microglia-PVM"   "L6b"             "Lamp5"          
# [17] "Sncg"            "L6 CT"           "L6 IT Car3"      "Chandelier"     
# [21] "Sst Chodl"       "Pax6"            "L5 ET"           "Endothelial"    

# Unique Liu MTC subclasses:
#  [1] "OPC"              "Oli"              "Ast"              "Exc L2-3 IT"     
#  [5] "Exc L4-5 IT-1"    "Inh SST"          "Exc L4-5 IT-2"    "Inh PAX6"        
#  [9] "Inh PVALB"        "Exc L3-4 IT"      "Exc L6b"          "SMC"             
# [13] "Mic"              "Exc L5/6 NP"      "Inh VIP"          "Exc L5/6 IT Car3"
# [17] "End"              "Exc L5-6 IT"      "Exc L6 IT"        "Exc L5 ET"       
# [21] "Inh LAMP5"        "Exc L3-5 IT"      "Exc L6 CT"        "Fib"             
# [25] "Per"              "Exc EC"           "T"                "Exc HC"          
# [29] "Epd"   