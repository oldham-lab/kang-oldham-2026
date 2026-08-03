# scz_db: Schizophrenia gene database for DFC brain region
# Rebuilt 2026-07-22 from dfc_overlaps.csv (HGNC-corrected) via pre-fetched-PMID Sonnet-5
#   literature review (deterministic NCBI prefetch + alias augmentation) + curated prefilter DB.
# n = number of SCZ-linked papers (up to 5; '5+' if more); n = 1 for curated-database genes
# ref = PMID of most relevant reference (literature genes); ref = website for database-sourced genes

scz_db_dfc <- list(

  # Synaptic vesicle cycling / neurotransmitter release
  RIMS2 = list(n = '1', ref = 18490030),  # RIMS2, a cytomatrix active-zone protein regulating synaptic vesicle release, shows incre...

  # Dendritic spine morphology or density
  ROCK2 = list(n = '3', ref = 41168987),  # ROCK2 hyperactivation drives dendritic spine loss and synaptic ultrastructural defects i...

  # Neurodevelopmental processes
  ACTR2 = list(n = '1', ref = 24996170),  # Genetic analysis found an epistatic interaction between CYFIP1 (15q11.2 CNV locus) and W...
  PAFAH1B1 = list(n = '1', ref = 16510495),  # LIS1/PAFAH1B1, a DISC1-interacting protein important for neuronal migration, shows reduc...
  SBNO1 = list(n = '2', ref = 36057539),  # SBNO1 is expressed in the developing cortical plate and required for neurite/axon-dendri...

  # Epigenetic regulation
  NSD2 = list(n = '1', ref = 32169559),  # Genetic interaction between histone-modification genes H3F3B and NSD2 (an H3K36 methyltr...
  PPM1E = list(n = '1', ref = 22832356),  # Promoter histone H3K9K14 acetylation at PPM1E is significantly reduced (hypoacetylated) ...

  # RNA processing / splicing / post-transcriptional regulation
  IPO5 = list(n = '1', ref = 20542336),  # Combined family/case-control genetic study found weak association of IPO5 with schizophr...
  UBE2K = list(n = '3', ref = 30901725),  # UBE2K shows elevated blood and brain protein/mRNA levels correlating with positive psych...

  # Schizophrenia GWAS (common-variant)
  NLGN1 = list(n = '3', ref = 26674772),  # Case-control genetic association study found NLGN1 polymorphisms (cell-adhesion-molecule...
  NTRK3 = list(n = '4', ref = 19344762),  # Multiple studies associate common NTRK3 (neurotrophin receptor) variants with schizophre...
  PCDH7 = list(n = '3', ref = 29503163),  # PCDH7 SNPs show genome-wide-significant association with antipsychotic treatment respons...

  # Calcium-channel signaling
  ATP2B1 = list(n = '1', ref = 42436150),  # Cross-trait GWAS pleiotropy analysis identified ATP2B1 as a schizophrenia risk locus inv...

  # Rare coding / CNV variant burden
  AKAP9 = list(n = '2', ref = 25943950),  # Case-only resequencing identified a rare damaging missense variant (K873R) in AKAP9 asso...
  ARID4A = list(n = '1', ref = 35365808),  # Whole-exome sequencing in Han Chinese schizophrenia-affected sibling families identified...
  CLTC = list(n = '1', ref = 36645932),  # Diagnostic clinical exome sequencing of childhood-onset schizophrenia families identifie...
  DNM3 = list(n = '1', ref = 26136833),  # Postmortem brain copy-number analysis in schizophrenia patients identified a gene-dosage...
  KPNA1 = list(n = '3', ref = 34767585),  # A de novo KPNA1 (importin alpha5) mutation has been linked to schizophrenia; KPNA1 regul...

  # Prefilter database entries
  AKAP6 = list(n = 1, ref = "platform.opentargets.org"),  # A-kinase anchoring protein 6; SCZ-associated in OpenTargets_Platform
  ARFGEF2 = list(n = 1, ref = "platform.opentargets.org"),  # ARF guanine nucleotide exchange factor 2; SCZ-associated in OpenTargets_Platform
  ATRX = list(n = 1, ref = "ncbi.nlm.nih.gov/clinvar"),  # ATRX chromatin remodeler; SCZ-associated in OpenTargets_Platform|ClinVar
  CACNA2D1 = list(n = 1, ref = "platform.opentargets.org"),  # calcium voltage-gated channel auxiliary subunit alpha2delta 1; SCZ-associated in OpenTargets_Platform
  CAMSAP2 = list(n = 1, ref = "platform.opentargets.org"),  # calmodulin regulated spectrin associated protein family member 2; SCZ-associated in OpenTargets_Platform
  CLCN4 = list(n = 1, ref = "ncbi.nlm.nih.gov/clinvar"),  # Cl-/H+ antiporter 4; SCZ-associated in OpenTargets_Platform|ClinVar
  CNOT1 = list(n = 1, ref = "platform.opentargets.org"),  # CCR4-NOT transcription complex subunit 1; SCZ-associated in OpenTargets_Platform
  CNTNAP2 = list(n = 1, ref = "ncbi.nlm.nih.gov/clinvar"),  # contactin associated protein 2; SCZ-associated in OpenTargets_Platform|ClinVar
  DENND5B = list(n = 1, ref = "platform.opentargets.org"),  # DENN domain containing 5B; SCZ-associated in OpenTargets_Platform
  DLG2 = list(n = 1, ref = "ncbi.nlm.nih.gov/clinvar"),  # discs large MAGUK scaffold protein 2; SCZ-associated in ClinVar
  EXOC4 = list(n = 1, ref = "platform.opentargets.org"),  # exocyst complex component 4; SCZ-associated in OpenTargets_Platform
  KIF21A = list(n = 1, ref = "platform.opentargets.org"),  # kinesin family member 21A; SCZ-associated in OpenTargets_Platform
  KIF5C = list(n = 1, ref = "platform.opentargets.org"),  # kinesin family member 5C; SCZ-associated in OpenTargets_Platform
  NALCN = list(n = 1, ref = "platform.opentargets.org"),  # sodium leak channel, non-selective; SCZ-associated in OpenTargets_Platform
  NBAS = list(n = 1, ref = "ncbi.nlm.nih.gov/clinvar"),  # NBAS subunit of NRZ tethering complex; SCZ-associated in OpenTargets_Platform|ClinVar
  NCOA7 = list(n = 1, ref = "platform.opentargets.org"),  # nuclear receptor coactivator 7; SCZ-associated in OpenTargets_Platform
  NETO1 = list(n = 1, ref = "platform.opentargets.org"),  # neuropilin and tolloid like 1; SCZ-associated in OpenTargets_Platform
  NSF = list(n = 1, ref = "omim.org"),  # N-ethylmaleimide sensitive factor, vesicle fusing ATPase; SCZ-associated in OMIM
  PBRM1 = list(n = 1, ref = "platform.opentargets.org"),  # polybromo 1; SCZ-associated in OpenTargets_Platform
  PIK3R1 = list(n = 1, ref = "platform.opentargets.org"),  # phosphoinositide-3-kinase regulatory subunit 1; SCZ-associated in OpenTargets_Platform
  PLXNA4 = list(n = 1, ref = "platform.opentargets.org"),  # plexin A4; SCZ-associated in OpenTargets_Platform
  PPIG = list(n = 1, ref = "platform.opentargets.org"),  # peptidylprolyl isomerase G; SCZ-associated in OpenTargets_Platform
  RAPGEF4 = list(n = 1, ref = "platform.opentargets.org"),  # Rap guanine nucleotide exchange factor 4; SCZ-associated in OpenTargets_Platform
  RB1CC1 = list(n = 1, ref = "ncbi.nlm.nih.gov/clinvar"),  # RB1 inducible coiled-coil 1; SCZ-associated in ClinVar
  RPS13 = list(n = 1, ref = "platform.opentargets.org"),  # ribosomal protein S13; SCZ-associated in OpenTargets_Platform
  SCN2A = list(n = 1, ref = "omim.org"),  # sodium voltage-gated channel alpha subunit 2; SCZ-associated in OMIM|OpenTargets_Platform
  SRPK2 = list(n = 1, ref = "platform.opentargets.org"),  # SRSF protein kinase 2; SCZ-associated in OpenTargets_Platform
  SYNE1 = list(n = 1, ref = "platform.opentargets.org"),  # spectrin repeat containing nuclear envelope protein 1; SCZ-associated in OpenTargets_Platform
  TCF4 = list(n = 1, ref = "platform.opentargets.org"),  # transcription factor 4; SCZ-associated in OpenTargets_Platform
  TGFBRAP1 = list(n = 1, ref = "platform.opentargets.org"),  # transforming growth factor beta receptor associated protein 1; SCZ-associated in OpenTargets_Platform
  TNKS = list(n = 1, ref = "platform.opentargets.org"),  # tankyrase; SCZ-associated in OpenTargets_Platform
  WDR7 = list(n = 1, ref = "platform.opentargets.org"),  # WD repeat domain 7; SCZ-associated in OpenTargets_Platform
  YTHDC2 = list(n = 1, ref = "platform.opentargets.org")  # YTH N6-methyladenosine RNA binding protein C2; SCZ-associated in OpenTargets_Platform
)
