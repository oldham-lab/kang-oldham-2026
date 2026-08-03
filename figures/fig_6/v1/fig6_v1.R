library(tidyverse)
library(data.table)
library(qs)

# Fig 6
# - dCOPA schematic
# - Examples of dCOPA modules

#####################
# Panel A - schematic
#####################

save_dir <- file.path(Sys.getenv("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_6/")

# # Real mod (con)
# real_mod_con <- data.frame("ct" = c("A", "B", "C"),
#                            "val" = c(1.2, 4.1, 2.5))
# p <- ggplot(real_mod_con, aes(x = ct, y = val)) +
#   theme_void() +
#   labs(x = "Celltype", y = "REI", title = "Control") +
#   geom_bar(stat = "identity", fill = "darkgrey") +
#   theme(plot.title = element_text(size = 18, hjust = 0.5),
#         axis.text.x = element_text(size = 18),
#         axis.title.x = element_text(size = 18),
#         axis.title.y = element_text(size = 18, angle = 90))
# ggsave(p, file = file.path(save_dir, "schematic", "real_mod_con.pdf"), width = 3, height = 2)

# # Real mod (disease)
# real_mod_dis <- data.frame("ct" = c("A", "B", "C", "A", "B", "C"),
#                            "type" = c("a", "a", "a", "b", "b", "b"),
#                            "val" = c(1.2, 3.2, 2.5, 0, .9, 0)) |>
#   mutate(type = factor(type, levels = rev(unique(type))))
# p <- ggplot(real_mod_dis, aes(x = ct, y = val, fill = type)) +
#   theme_void() +
#   labs(x = "Celltype", y = "REI", title = "Disease") +
#   geom_bar(stat = "identity", position = "stack") +
#   theme(plot.title = element_text(size = 18, hjust = 0.5),
#         axis.text.x = element_text(size = 18),
#         axis.title.x = element_text(size = 18),
#         axis.title.y = element_text(size = 18, angle = 90),
#         legend.position = "none") +
#   scale_fill_manual(values = c("a" = "darkgrey", "b" = "red")) 
# ggsave(p, file = file.path(save_dir, "schematic", "real_mod_dis.pdf"), width = 3, height = 2)

# Real mod
real_mod <- data.frame("ct" = c("A", "B", "C", "A", "B", "C"),
                       "type" = c("Control", "Control", "Control", "Disease", "Disease", "Disease"),
                       "val" = c(1.2, 2.5, 2.5, 1.2, 4.1, 2.5),
                       "lab" = c("", "*", "", "", "*", "")) 
p <- ggplot(real_mod, aes(x = ct, y = val, fill = type, label = lab)) +
  theme_void() +
  labs(x = "Celltype", y = "REI") +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(size = 8, vjust = -0.5) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18, angle = 90),
        legend.position = "none") +
  scale_fill_manual(values = c("Control" = "#EAEBEB", "Disease" = "red")) 
ggsave(p, file = file.path(save_dir, "schematic", "real_mod.pdf"), width = 3, height = 2)


# Density
dmat <- data.frame("vals" = rnorm(1000))
p <- ggplot(dmat, aes(x = vals)) + 
  theme_void() +
  geom_density(fill = "grey") +
  geom_vline(xintercept = 2.67, color = "red", linetype = "dashed")
ggsave(p, file = file.path(save_dir, "schematic", "density.pdf"), width = 3, height = 2)

