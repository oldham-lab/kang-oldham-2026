# scz_db: Schizophrenia gene database for DFC brain region
# Built from dfc_overlaps.csv via three-agent literature pipeline
# n = number of SCZ-linked papers (up to 5; "5+" if more);
#     n = 1 for genes sourced from a curated SCZ database (no literature review run)
# ref = PMID of most relevant verified reference (literature-pipeline genes);
#     ref = "website" for database-sourced genes (e.g., "omim.org", "platform.opentargets.org")

scz_db_dfc <- list(

  # Dopaminergic signaling dysfunction
  GPR12 = list(n = '3', ref = 37605621),  # GPR12 is a brain-specific orphan GPCR identified as an emerging therap

  # Synaptic vesicle cycling or neurotransmitter release
  DNM3 = list(n = '2', ref = 26136833),  # Dynamin 3 (DNM3) shows copy number variation in brain tissue of schizo
  ENTPD4 = list(n = '1', ref = 20537721),  # ENTPD4 polymorphisms were directly tested for schizophrenia susceptibi
  RIMS2 = list(n = '2', ref = 18490030),  # RIMS2 expression is significantly increased in the amygdala of schizop
  UBE2K = list(n = '3', ref = 30901725),  # UBE2K protein levels are significantly elevated in blood and brain of
  UBE2N = list(n = '1', ref = 12363385),  # UBE2N (ubiquitin-conjugating enzyme E2N) shows altered expression in t
  UNC13C = list(n = '2', ref = 38188011),  # UNC13C (Munc13-3) is a presynaptic priming factor for synaptic vesicle

  # Dendritic spine morphology or density
  ROCK2 = list(n = '3', ref = 41922796),  # ROCK2 regulates actin cytoskeleton dynamics and dendritic spine morpho

  # Neurodevelopmental processes
  AJAP1 = list(n = '1', ref = 19850283),  # Novel candidate gene for schizophrenia susceptibility and risperidone
  AKAP9 = list(n = '3', ref = 25943950),  # AKAP9 harbors rare coding variants significantly associated with schiz
  CNTNAP5 = list(n = '5', ref = 29503163),  # CNTNAP5 was identified as one of five novel loci in a GWAS of antipsyc
  KPNA1 = list(n = '4', ref = 38336912),  # De novo mutations in KPNA1 (importin α5) are linked to schizophrenia;
  LRRTM3 = list(n = '1', ref = 24587117),  # LRRTM3 belongs to the LRRTM family of postsynaptic adhesion proteins t
  NLGN1 = list(n = '5+', ref = 26674772),  # NLGN1 polymorphisms are directly associated with schizophrenia suscept
  NTRK3 = list(n = '5+', ref = 19344762),  # NTRK3 (TrkC) polymorphisms are associated with hippocampal function an
  PAFAH1B1 = list(n = '4', ref = 16510495),  # PAFAH1B1 (LIS1), a critical regulator of neuronal migration and dynein
  PCDH7 = list(n = '5+', ref = 29503163),  # PCDH7 was identified as a novel locus in a GWAS of antipsychotic treat
  SBNO1 = list(n = '3', ref = 36057539),  # SBNO1 (strawberry notch homolog 1) regulates neurite outgrowth in cort

  # Mitochondrial dysfunction
  PGRMC1 = list(n = '5+', ref = 37087062),  # PGRMC1 (progesterone receptor membrane component 1 / sigma-2 receptor)

  # Epigenetic regulation
  ARID4A = list(n = '1', ref = 35365808),  # A novel heterozygous missense variant in ARID4A was identified in Han
  NSD2 = list(n = '1', ref = 32169559),  # NSD2 encodes a histone H3K36 dimethyltransferase whose variants intera

  # RNA processing or post-transcriptional regulation
  IPO5 = list(n = '1', ref = 20542336),  # IPO5 (importin 5) shows genetic association with schizophrenia in fami

  # Prefilter database entries
  # Source: OMIM
  NSF = list(n = 1, ref = "omim.org"),  # N-ethylmaleimide sensitive factor, vesicle fusing ATPase; SCZ-associated in OMIM
  SCN2A = list(n = 1, ref = "omim.org"),  # sodium voltage-gated channel alpha subunit 2; SCZ-associated in OMIM; OpenTargets_Platform

  # Source: ClinVar
  ATRX = list(n = 1, ref = "ncbi.nlm.nih.gov/clinvar"),  # ATRX chromatin remodeler; SCZ-associated in OpenTargets_Platform; ClinVar
  CLCN4 = list(n = 1, ref = "ncbi.nlm.nih.gov/clinvar"),  # Cl-/H+ antiporter 4; SCZ-associated in OpenTargets_Platform; ClinVar
  CNTNAP2 = list(n = 1, ref = "ncbi.nlm.nih.gov/clinvar"),  # contactin associated protein 2; SCZ-associated in OpenTargets_Platform; ClinVar
  DLG2 = list(n = 1, ref = "ncbi.nlm.nih.gov/clinvar"),  # discs large MAGUK scaffold protein 2; SCZ-associated in ClinVar
  NBAS = list(n = 1, ref = "ncbi.nlm.nih.gov/clinvar"),  # NBAS subunit of NRZ tethering complex; SCZ-associated in OpenTargets_Platform; ClinVar
  RB1CC1 = list(n = 1, ref = "ncbi.nlm.nih.gov/clinvar"),  # RB1 inducible coiled-coil 1; SCZ-associated in ClinVar

  # Source: OpenTargets_Platform
  AKAP6 = list(n = 1, ref = "platform.opentargets.org"),  # A-kinase anchoring protein 6; SCZ-associated in OpenTargets_Platform
  ARFGEF2 = list(n = 1, ref = "platform.opentargets.org"),  # ARF guanine nucleotide exchange factor 2; SCZ-associated in OpenTargets_Platform
  CACNA2D1 = list(n = 1, ref = "platform.opentargets.org"),  # calcium voltage-gated channel auxiliary subunit alpha2delta 1; SCZ-associated in OpenTargets_Platform
  CAMSAP2 = list(n = 1, ref = "platform.opentargets.org"),  # calmodulin regulated spectrin associated protein family member 2; SCZ-associated in OpenTargets_Platform
  CNOT1 = list(n = 1, ref = "platform.opentargets.org"),  # CCR4-NOT transcription complex subunit 1; SCZ-associated in OpenTargets_Platform
  COQ10B = list(n = 1, ref = "platform.opentargets.org"),  # coenzyme Q10B; SCZ-associated in OpenTargets_Platform
  DENND5B = list(n = 1, ref = "platform.opentargets.org"),  # DENN domain containing 5B; SCZ-associated in OpenTargets_Platform
  EXOC4 = list(n = 1, ref = "platform.opentargets.org"),  # exocyst complex component 4; SCZ-associated in OpenTargets_Platform
  KIF21A = list(n = 1, ref = "platform.opentargets.org"),  # kinesin family member 21A; SCZ-associated in OpenTargets_Platform
  KIF5C = list(n = 1, ref = "platform.opentargets.org"),  # kinesin family member 5C; SCZ-associated in OpenTargets_Platform
  NALCN = list(n = 1, ref = "platform.opentargets.org"),  # sodium leak channel, non-selective; SCZ-associated in OpenTargets_Platform
  NCOA7 = list(n = 1, ref = "platform.opentargets.org"),  # nuclear receptor coactivator 7; SCZ-associated in OpenTargets_Platform
  NETO1 = list(n = 1, ref = "platform.opentargets.org"),  # neuropilin and tolloid like 1; SCZ-associated in OpenTargets_Platform
  PBRM1 = list(n = 1, ref = "platform.opentargets.org"),  # polybromo 1; SCZ-associated in OpenTargets_Platform
  PIK3R1 = list(n = 1, ref = "platform.opentargets.org"),  # phosphoinositide-3-kinase regulatory subunit 1; SCZ-associated in OpenTargets_Platform
  PLXNA4 = list(n = 1, ref = "platform.opentargets.org"),  # plexin A4; SCZ-associated in OpenTargets_Platform
  PPIG = list(n = 1, ref = "platform.opentargets.org"),  # peptidylprolyl isomerase G; SCZ-associated in OpenTargets_Platform
  RAPGEF4 = list(n = 1, ref = "platform.opentargets.org"),  # Rap guanine nucleotide exchange factor 4; SCZ-associated in OpenTargets_Platform
  RPS13 = list(n = 1, ref = "platform.opentargets.org"),  # ribosomal protein S13; SCZ-associated in OpenTargets_Platform
  SRPK2 = list(n = 1, ref = "platform.opentargets.org"),  # SRSF protein kinase 2; SCZ-associated in OpenTargets_Platform
  SYNE1 = list(n = 1, ref = "platform.opentargets.org"),  # spectrin repeat containing nuclear envelope protein 1; SCZ-associated in OpenTargets_Platform
  TCF4 = list(n = 1, ref = "platform.opentargets.org"),  # transcription factor 4; SCZ-associated in OpenTargets_Platform
  TGFBRAP1 = list(n = 1, ref = "platform.opentargets.org"),  # transforming growth factor beta receptor associated protein 1; SCZ-associated in OpenTargets_Platform
  TNKS = list(n = 1, ref = "platform.opentargets.org"),  # tankyrase; SCZ-associated in OpenTargets_Platform
  WDR7 = list(n = 1, ref = "platform.opentargets.org"),  # WD repeat domain 7; SCZ-associated in OpenTargets_Platform
  YTHDC2 = list(n = 1, ref = "platform.opentargets.org")  # YTH N6-methyladenosine RNA binding protein C2; SCZ-associated in OpenTargets_Platform

)
