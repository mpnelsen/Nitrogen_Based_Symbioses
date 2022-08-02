# Nitrogen_Based_Symbioses
Data and trees used for stochastic mapping, and timescale used for plotting in Boyce et al. (2022) Geobiology

Files from previous work were retrieved from the following sources:

Agaricomycotina phylogenies of Varga et al. (2019 Nature Ecology & Evolution) were downloaded from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.gc2k9r9

Trait data of Tedersoo et al. (2020 Frontiers in Microbiology) used for Agaricomycotina analyses were downloaded from: https://plutof.ut.ee/#/doi/10.15156/BIO/807446

Lecanoromycetes data matrix from Nelsen et al. (2020 PNAS) downloaded from: https://github.com/mpnelsen/Lecanoromycetes_megaphylogeny/blob/master/HiSSE/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_anycyano_forcorhmm.csv

Lecanoromycetes topologies from Nelsen et al. (2020 PNAS) downloaded from: https://github.com/mpnelsen/Lecanoromycetes_megaphylogeny/tree/master/Trees

International Chronostratigraphic Chart (ICS 2020) was transcribed from: https://stratigraphy.org/ICSchart/ChronostratChart2020-03.pdf

Files included here:

R script used for analyses:
code_24_clean2.r

Agaricomycotina species coded for whether or not parent genus forms ectomycorrhizal associations in Tedersoo et al. (2020):
EctoStates_File.csv

Agaricomycotina phylogenies from Varga et al. (2019) with names modified, and reduced to unique, formally-described species in genera that were clearly ectomycorrhizal or not:
Ecto_Trees_Cleaned.tre

Lecanoromycetes species coded for whether or not they associate with cyanobacteria from Nelsen et al. (2020):
Lecanoromycetes_anycyano.csv

Lecanoromycetes phylogenies used for analyses from Nelsen et al. (2020).  Trees 1-10 are derived from bootstrap replicates, and tree 11 is the ML tree:
Lecanoromycetes_10boots_1ML.tre

ICS 2020 timescale used for plotting
timescale_ics2020.csv

Changes and change rates in 10my timebins for all trees analyzed:
lec_change_10my.csv
agarico_change_10my.csv

Changes and change rates in 10my timebins for Lecanoromycetes ML tree and single Agaricomycotina tree:
lec_change_10my_ML_ONLY.csv
agarico_change_10my_Tree4_53_ONLY.csv

