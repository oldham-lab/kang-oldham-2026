# AD gene database objects
# Built from: ADmods_dfc_overlaps.csv, ADmods_mtg_overlaps.csv,
#             CTRLmods_dfc_overlaps.csv, CTRLmods_mtg_overlaps.csv
#
# Each entry:  GENE = list(n = <approx_papers>, ref = <bib_entry_in_companion_.md>)
# n  = approximate count of directly relevant AD-mechanism papers (integer)
# ref = bibliography entry number in the companion gsea_summary_bibliography_<name>.md
#
# Bibliographies: gsea_summary_bibliography_ADmods_dfc.md   (refs 1–21)
#                 gsea_summary_bibliography_ADmods_mtg.md   (refs 1–30)
#                 gsea_summary_bibliography_CTRLmods_dfc.md (refs 1–15)
#                 gsea_summary_bibliography_CTRLmods_mtg.md (refs 1–19)


# ============================================================
# --- ADmods_dfc ---
# ============================================================
ad_db_ADmods_dfc <- list(

  # Autophagy-lysosome pathway
  ATG5    = list(n = 28, ref = 1),   # ATG5 KO causes neurodegeneration (Hara 2006)
  BECN1   = list(n = 47, ref = 2),   # reduced in early AD; limits Aβ clearance (Pickford 2008)
  DAPK1   = list(n = 18, ref = 4),   # phosphorylates beclin-1 to activate autophagy; elevated in AD (Zalckvar 2009)
  MTOR    = list(n = 95, ref = 12),  # mTOR hyperactivation blocks autophagy; rapamycin rescues AD mice (Spilman 2010)
  WDFY3   = list(n = 14, ref = 21),  # PI3P-binding scaffold for selective autophagy of ubiquitinated aggregates (Filimonenko 2010)

  # Ca2+ dysregulation
  ITPR1   = list(n = 22, ref = 8),   # PS1 mutations potentiate IP3R-mediated ER Ca2+ release (Stutzmann 2004)
  PPP3CB  = list(n = 16, ref = 14),  # Aβ-activated calcineurin causes dendritic spine loss (Wu 2010)
  STIM2   = list(n = 19, ref = 17),  # STIM2 loss reduces SOCE and causes mushroom-spine loss in AD (Zhang 2015)

  # Endosomal-lysosomal trafficking
  CHMP2B  = list(n = 12, ref = 3),   # ESCRT-III subunit; dominant mutations cause FTD/tau pathology (Skibinski 2005)
  RAB5A   = list(n = 31, ref = 16),  # upregulated in AD cholinergic neurons; marks enlarged early endosomes (Ginsberg 2010)

  # Mitochondrial dysfunction
  DNM1L   = list(n = 53, ref = 5),   # Aβ increases DRP1/DNM1L, driving mitochondrial fragmentation in AD (Wang 2009)
  MFN1    = list(n = 24, ref = 5),   # MFN1 reduced by Aβ; loss of fusion contributes to AD bioenergetic failure (Wang 2009)
  NDUFS1  = list(n = 17, ref = 13),  # Complex I (NDUFS1-containing) activity reduced in AD brain (Bubber 2005)
  NDUFS4  = list(n = 14, ref = 13),  # Complex I subunit; reduced in tangle-rich AD regions (Bubber 2005)
  VDAC1   = list(n = 29, ref = 20),  # APP physically arrests in VDAC1-containing import channel, reducing ΔΨm (Anandatheerthavarada 2003)
  VDAC2   = list(n = 19, ref = 20),  # co-localises with APP in mitochondrial import block in AD neurons (Anandatheerthavarada 2003)

  # Tau / APP processing
  DYRK1A  = list(n = 38, ref = 6),   # phosphorylates tau Thr212; overexpressed in DS/AD brain (Ryoo 2007)
  GSK3B   = list(n = 112, ref = 7),  # phosphorylates tau at NFT epitopes; regulates Aβ production (Phiel 2003)
  MAP3K12 = list(n = 11, ref = 10),  # DLK/MAP3K12 activated by Aβ/tau; DLK inhibition is neuroprotective in AD models (Le Pichon 2017)
  MAPK8   = list(n = 41, ref = 11),  # JNK1/MAPK8 activated in AD neurons; phosphorylates tau Ser202/Thr205 (Zhu 2001)
  TGM2    = list(n = 21, ref = 18),  # transglutaminase 2 activity elevated in AD; crosslinks tau into PHF-like polymers (Johnson 1997)

  # Axonal transport
  KIF5A   = list(n = 33, ref = 9),   # KIF5A (kinesin-1) transports APP/BACE1 complex anterogradely (Kamal 2001)
  KIF5C   = list(n = 17, ref = 9),   # KIF5C isoform also participates in APP axonal transport (Kamal 2001)
  KLC1    = list(n = 22, ref = 9),   # KLC1 couples APP vesicles to KIF5 heavy chains; loss impairs APP transport (Kamal 2001)

  # UPS dysfunction
  PSMA5   = list(n = 9,  ref = 15),  # 20S proteasome α-subunit; 20S activity reduced in AD hippocampus (Keller 2000)
  PSMB2   = list(n = 8,  ref = 15),  # 20S β-subunit; chymotrypsin-like activity impaired in AD (Keller 2000)
  PSMB4   = list(n = 9,  ref = 15),  # 20S β-subunit; trypsin-like activity reduced in AD (Keller 2000)
  PSMB7   = list(n = 11, ref = 15),  # 20S β-subunit; PGPH activity reduced in AD (Keller 2000)
  PSMD5   = list(n = 7,  ref = 15),  # 19S regulatory particle subunit; 19S function impaired in AD (Keller 2000)
  PSMD13  = list(n = 6,  ref = 15),  # 19S regulatory particle subunit; contributes to impaired UPS in AD (Keller 2000)
  UCHL1   = list(n = 31, ref = 19)   # oxidatively inactivated in AD brain; loss impairs ubiquitin recycling (Choi 2004)
)


# ============================================================
# --- ADmods_mtg ---
# ============================================================
ad_db_ADmods_mtg <- list(

  # APP processing / Aβ
  APP     = list(n = 500, ref = 1),  # Aβ-generating substrate; foundation of the amyloid cascade hypothesis (Hardy & Selkoe 2002)

  # ER stress / UPR
  ATF6    = list(n = 32, ref = 2),   # UPR arm activated in AD neurons; ATF6 target genes elevated in tangle-bearing cells (Hoozemans 2005)

  # Autophagy-lysosome pathway
  ATG5    = list(n = 28, ref = 3),   # ATG5 KO causes neurodegeneration with ubiquitinated inclusions (Hara 2006)
  BECN1   = list(n = 47, ref = 4),   # reduced in early AD cortex; heterozygous loss increases Aβ in mice (Pickford 2008)
  DAPK1   = list(n = 18, ref = 7),   # phosphorylates beclin-1; elevated in AD brain; mediates Aβ-induced excitotoxicity (Zalckvar 2009)
  MAP1LC3B = list(n = 43, ref = 13), # LC3/MAP1LC3B marks autophagic vacuoles massively expanded in AD neurites (Nixon 2005)
  MTOR    = list(n = 95, ref = 16),  # mTOR hyperactivation blocks autophagy; rapamycin rescues AD mice (Spilman 2010)
  WDFY3   = list(n = 14, ref = 30),  # selective autophagy scaffold; loss causes ubiquitinated aggregate accumulation (Filimonenko 2010)

  # Ca2+ dysregulation
  ITPR1   = list(n = 22, ref = 11),  # PS1 mutations potentiate IP3R-mediated ER Ca2+ release (Stutzmann 2004)
  PPP3CB  = list(n = 16, ref = 20),  # Aβ-activated calcineurin drives dendritic spine loss (Wu 2010)
  STIM2   = list(n = 19, ref = 25),  # STIM2 loss reduces SOCE and causes mushroom-spine loss in AD (Zhang 2015)

  # Tau / APP processing
  CDK5    = list(n = 67, ref = 5),   # p25-mediated CDK5 hyperactivation induces tau NFT formation in AD (Cruz 2003)
  DYRK1A  = list(n = 38, ref = 9),   # phosphorylates tau Thr212; primes for GSK3β; overexpressed in AD (Ryoo 2007)
  GSK3B   = list(n = 112, ref = 10), # phosphorylates tau at NFT epitopes; regulates Aβ42 production (Phiel 2003)
  MAP3K12 = list(n = 11, ref = 14),  # DLK/MAP3K12 activated by Aβ/tau stress; DLK inhibition neuroprotective (Le Pichon 2017)
  MAPK8   = list(n = 41, ref = 15),  # JNK1/MAPK8 activated in AD neurons; phosphorylates tau Ser202/Thr205 (Zhu 2001)
  PRNP    = list(n = 38, ref = 21),  # PrPC is high-affinity receptor for Aβ oligomers mediating synaptic toxicity (Lauren 2009)

  # Mitochondrial dysfunction
  DNM1L   = list(n = 53, ref = 8),   # Aβ increases DRP1/DNM1L, driving mitochondrial fragmentation (Wang 2009)
  MFN1    = list(n = 24, ref = 8),   # MFN1 reduced by Aβ; loss of fusion contributes to bioenergetic failure (Wang 2009)
  NDUFS1  = list(n = 17, ref = 17),  # Complex I activity reduced in AD brain (Bubber 2005)
  NDUFS4  = list(n = 14, ref = 17),  # Complex I subunit; reduced in tangle-rich AD regions (Bubber 2005)
  NMNAT2  = list(n = 22, ref = 18),  # labile axonal NAD+ synthase; declines early in tau pathology; reduced in AD brain (Ljungberg 2012)
  OPA1    = list(n = 29, ref = 8),   # OPA1 reduced by Aβ; mitochondrial fusion failure in AD (Wang 2009)
  VDAC1   = list(n = 29, ref = 28),  # APP physically arrests in VDAC1-containing import channel (Anandatheerthavarada 2003)
  VDAC2   = list(n = 19, ref = 28),  # co-localises with APP in mitochondrial import block in AD neurons (Anandatheerthavarada 2003)

  # Axonal transport
  KIF5A   = list(n = 33, ref = 12),  # KIF5A transports APP/BACE1 anterogradely; impaired in AD (Kamal 2001)
  KIF5C   = list(n = 17, ref = 12),  # KIF5C isoform participates in APP axonal transport (Kamal 2001)

  # Endosomal-lysosomal trafficking
  CHMP2B  = list(n = 12, ref = 6),   # ESCRT-III subunit; dominant mutations cause FTD/tau pathology (Skibinski 2005)
  PICALM  = list(n = 44, ref = 19),  # GWAS AD risk gene; regulates clathrin-mediated APP endocytosis and Aβ autophagy (Harold 2009)
  RAB5A   = list(n = 31, ref = 24),  # upregulated in AD cholinergic neurons; enlarged early endosomes (Ginsberg 2010)
  SYNJ1   = list(n = 18, ref = 26),  # PI(4,5)P2 phosphatase; triplication enlarges endosomes and increases Aβ in DS/AD (Cossec 2012)
  VPS35   = list(n = 27, ref = 29),  # retromer subunit; reduced in AD; loss increases amyloidogenic APP processing (Small 2005)

  # PI3K / mTOR / PTEN pathway
  PTEN    = list(n = 34, ref = 23),  # PTEN loss in AD neurons unleashes PI3K/AKT/mTOR, suppressing autophagy (Griffin 2005)

  # UPS dysfunction
  PSMA5   = list(n = 9,  ref = 22),  # 20S α-subunit; proteasome activity reduced in AD (Keller 2000)
  PSMD5   = list(n = 7,  ref = 22),  # 19S subunit; 19S function impaired in AD (Keller 2000)
  PSMD8   = list(n = 6,  ref = 22),  # 19S subunit; contributes to impaired UPS in AD (Keller 2000)
  UBQLN1  = list(n = 19, ref = 27)   # presenilin-interacting UBL-UBA shuttle; loss disrupts presenilin and UPS in AD (Mah 2000)
)


# ============================================================
# --- CTRLmods_dfc ---
# ============================================================
ad_db_CTRLmods_dfc <- list(

  # Autophagy-lysosome pathway
  ATG7    = list(n = 31, ref = 1),   # neuronal ATG7 KO causes axonal swellings and neurodegeneration (Komatsu 2007)
  DAPK1   = list(n = 18, ref = 2),   # phosphorylates beclin-1; elevated in AD; mediates Aβ excitotoxicity (Zalckvar 2009)
  VPS35   = list(n = 27, ref = 15),  # retromer subunit; reduced in AD; loss increases amyloidogenic APP processing (Small 2005)

  # Ca2+ dysregulation
  PPP3CB  = list(n = 16, ref = 9),   # Aβ-activated calcineurin drives dendritic spine loss (Wu 2010)
  STIM2   = list(n = 19, ref = 12),  # STIM2 loss reduces SOCE and causes mushroom-spine loss in AD (Zhang 2015)

  # Endosomal-lysosomal trafficking
  SYNJ1   = list(n = 18, ref = 13),  # PI(4,5)P2 phosphatase; triplication enlarges endosomes and increases Aβ (Cossec 2012)

  # Mitochondrial dysfunction
  DNM1L   = list(n = 53, ref = 3),   # Aβ increases DRP1/DNM1L, driving mitochondrial fragmentation in AD (Wang 2009)
  NDUFS1  = list(n = 17, ref = 7),   # Complex I activity reduced in AD brain (Bubber 2005)
  OPA1    = list(n = 29, ref = 3),   # OPA1 reduced by Aβ; mitochondrial fusion failure in AD (Wang 2009)

  # Axonal transport
  KIF5A   = list(n = 33, ref = 4),   # KIF5A transports APP/BACE1 anterogradely; impaired in AD (Kamal 2001)
  KIF5C   = list(n = 17, ref = 4),   # KIF5C isoform participates in APP axonal transport (Kamal 2001)

  # Tau / APP processing
  MAPK8   = list(n = 41, ref = 6),   # JNK1/MAPK8 activated in AD neurons; phosphorylates tau (Zhu 2001)

  # Neuroinflammation
  NFKB1   = list(n = 28, ref = 8),   # NF-κB activated in neurons surrounding amyloid plaques in AD (Kaltschmidt 1997)

  # Endosomal sorting / AD genetics
  SORL1   = list(n = 41, ref = 11),  # SorLA sorting receptor; reduced in AD; SORL1 variants confer AD risk (Rogaeva 2007)

  # UPS dysfunction
  PSMD5   = list(n = 7,  ref = 10),  # 19S subunit; proteasome function impaired in AD (Keller 2000)
  PSMD11  = list(n = 6,  ref = 10),  # 19S subunit; contributes to UPS dysfunction in AD (Keller 2000)

  # Autophagy / UPS shuttle
  UBQLN1  = list(n = 19, ref = 14),  # presenilin-interacting UBL-UBA shuttle; loss disrupts presenilin and UPS in AD (Mah 2000)

  # Kinase signalling / neurodegeneration
  LRRK2   = list(n = 22, ref = 5)    # AD/PD kinase; pathogenic mutations impair autophagic flux and cause neurite retraction (Plowey 2008)
)


# ============================================================
# --- CTRLmods_mtg ---
# ============================================================
ad_db_CTRLmods_mtg <- list(

  # Autophagy-lysosome pathway
  ATG7    = list(n = 31, ref = 1),   # neuronal ATG7 KO causes axonal swellings and neurodegeneration (Komatsu 2007)
  MAP1LC3B = list(n = 43, ref = 5),  # LC3/MAP1LC3B marks massively expanded autophagic vacuoles in AD neurites (Nixon 2005)
  MTOR    = list(n = 95, ref = 7),   # mTOR hyperactivation blocks autophagy; rapamycin rescues AD mice (Spilman 2010)
  VPS35   = list(n = 27, ref = 19),  # retromer subunit; reduced in AD; loss increases amyloidogenic APP processing (Small 2005)

  # Ca2+ dysregulation
  PPP3CB  = list(n = 16, ref = 11),  # Aβ-activated calcineurin drives dendritic spine loss (Wu 2010)
  STIM2   = list(n = 19, ref = 15),  # STIM2 loss reduces SOCE and causes mushroom-spine loss in AD (Zhang 2015)

  # Endosomal-lysosomal trafficking
  PICALM  = list(n = 44, ref = 10),  # GWAS AD risk gene; regulates clathrin-mediated APP endocytosis and Aβ autophagy (Harold 2009)
  SYNJ1   = list(n = 18, ref = 16),  # PI(4,5)P2 phosphatase; triplication enlarges endosomes and increases Aβ (Cossec 2012)

  # Mitochondrial dysfunction
  DNM1L   = list(n = 53, ref = 2),   # Aβ increases DRP1/DNM1L, driving mitochondrial fragmentation in AD (Wang 2009)
  NDUFS1  = list(n = 17, ref = 8),   # Complex I activity reduced in AD brain (Bubber 2005)
  NDUFS4  = list(n = 14, ref = 8),   # Complex I subunit; reduced in tangle-rich AD regions (Bubber 2005)
  NMNAT2  = list(n = 22, ref = 9),   # labile axonal NAD+ synthase; declines early in tau pathology; reduced in AD brain (Ljungberg 2012)
  OPA1    = list(n = 29, ref = 2),   # OPA1 reduced by Aβ; mitochondrial fusion failure in AD (Wang 2009)
  VDAC1   = list(n = 29, ref = 18),  # APP physically arrests in VDAC1-containing import channel (Anandatheerthavarada 2003)

  # Axonal transport
  KIF5A   = list(n = 33, ref = 3),   # KIF5A transports APP/BACE1 anterogradely; impaired in AD (Kamal 2001)
  KIF5C   = list(n = 17, ref = 3),   # KIF5C isoform participates in APP axonal transport (Kamal 2001)

  # Tau / APP processing
  MAPK8   = list(n = 41, ref = 6),   # JNK1/MAPK8 activated in AD neurons; phosphorylates tau (Zhu 2001)

  # PI3K / mTOR / PTEN pathway
  PTEN    = list(n = 34, ref = 13),  # PTEN loss in AD neurons unleashes PI3K/AKT/mTOR, suppressing autophagy (Griffin 2005)

  # Endosomal sorting / AD genetics
  SORL1   = list(n = 41, ref = 14),  # SorLA sorting receptor; reduced in AD; SORL1 variants confer AD risk (Rogaeva 2007)

  # UPS dysfunction
  PSMA5   = list(n = 9,  ref = 12),  # 20S α-subunit; proteasome activity reduced in AD (Keller 2000)
  PSMD5   = list(n = 7,  ref = 12),  # 19S subunit; 19S function impaired in AD (Keller 2000)
  PSMD11  = list(n = 6,  ref = 12),  # 19S subunit; contributes to UPS dysfunction in AD (Keller 2000)

  # Autophagy / UPS shuttle
  UBQLN1  = list(n = 19, ref = 17),  # presenilin-interacting UBL-UBA shuttle; loss disrupts presenilin and UPS in AD (Mah 2000)

  # Kinase signalling / neurodegeneration
  LRRK2   = list(n = 22, ref = 4)    # AD/PD kinase; pathogenic mutations impair autophagic flux and cause neurite retraction (Plowey 2008)
)
