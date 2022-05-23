## code to prepare `bat_phylogeny` dataset goes here

library(ape)

mormop_string <- "
(((Ptda, Ptpe), Ptpa), Mome);
"
mormop_tree <- read.tree(text = mormop_string)
mormop_tree$edge.length <- c(17.06, 3.06, 14.82, 14.82, 17.88, 34.94)
# plot(mormop_tree)
# dist.nodes(mormop_tree)

phyllo_string <- "
(Mor, (Maca, (Dero, (Leye, (Arja, (Stli, Stlu))))));
"
phyllo_tree <- read.tree(text = phyllo_string)
phyllo_tree$edge.length <- c(4.9, 4.85, 34.99, 1.76, 33.23, 4.87, 28.36, 9.81, 18.55, 12.84, 5.71, 5.71)
# plot(phyllo_tree)
# dist.nodes(phyllo_tree)

mp_tree <- bind.tree(phyllo_tree, mormop_tree, where = 1)
# plot(mp_tree)
# dist.nodes(mp_tree)

emball_string <- "
(MP, Bapl);
"
emball_tree <- read.tree(text = emball_string)
emball_tree$edge.length <- c(14.82, 54.66)
emp_tree <- bind.tree(emball_tree, mp_tree, where = 1)
# plot(emp_tree)
# dist.nodes(emp_tree)

vesper_string <- "
(((Myyu, Myvo), Idph), ((Pihe, Epfu), (Anpa, (Laxa, (Laci, Labl)))));
"
vesper_tree <- read.tree(text = vesper_string)
vesper_tree$edge.length <- c(2.45, 20.2, 13.29, 13.29, 33.49, 3.73, 1.04, 31.17, 31.17, 0.97, 31.24, 11.45, 19.79, 1.96, 17.83, 17.83)
# plot(vesper_tree)
# dist.nodes(vesper_tree)

moloss_string <- "
(EMP, (Ves, (Nyfe, Tabr)));
"
moloss_tree <- read.tree(text = moloss_string)
moloss_tree$edge.length <- c(1.38, 6.01, 14.09, 30.1, 19.93, 19.93)
# plot(moloss_tree)
# dist.nodes(moloss_tree)

bat_phylogeny <- bind.tree(moloss_tree, vesper_tree, where = 2)
bat_phylogeny <- bind.tree(bat_phylogeny, emp_tree, where = 1)
plot(bat_phylogeny)
# dist.nodes(bat_phylogeny)

usethis::use_data(bat_phylogeny, overwrite = TRUE)
