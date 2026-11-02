## =============================================================
## 5_figures.R
## Figures for the log(CPD)-based manuscript, drawn at print size
## (6.5 in wide), 600 dpi. No grey: blue = significantly steeper than
## average / main series, red = significantly shallower, black = not
## significant or neutral marks. Run after 3_scaling_checks.R.
## =============================================================
library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
source("config.R")
FIG <- "figures"; dir.create(FIG, showWarnings = FALSE)
BLUE <- "#1F5FA8"; RED <- "#B8323A"; INK <- "#000000"; GREEN <- "#2E7D32"; PURPLE <- "#6A3D9A"
W <- 6.5; DPI <- 800
th <- theme_minimal(base_size = 10.5) +
  theme(panel.grid = element_blank(), axis.line = element_line(color = INK, linewidth = 0.45),
        axis.ticks = element_line(color = INK, linewidth = 0.45), axis.ticks.length = unit(3, "pt"),
        axis.text = element_text(color = INK, size = 10), axis.title = element_text(color = INK, size = 10.5),
        legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 10, color = INK),
        plot.title = element_text(face = "bold", size = 12, color = INK, margin = margin(0, 0, 6, 0)),
        plot.title.position = "plot", plot.margin = margin(4, 8, 4, 4))
dircol <- c(above = BLUE, below = RED, `n.s.` = INK)
dirlab <- c(above = "Steeper than average", below = "Shallower than average", `n.s.` = "Not significant")
inside_legend <- theme(legend.position = "bottom", legend.justification = "center", legend.box.just = "left",
                       legend.background = element_rect(fill = "white", color = INK, linewidth = 0.3),
                       legend.margin = margin(4, 6, 4, 6))
num <- jsonlite::fromJSON(O("lc_numbers.json"))

## ---- Figure 1: elasticity across specifications ----------------------
main <- fread(O("lc_table_main.csv")); rob <- fread(O("lc_table_robustness.csv"))
g <- function(m, t, c) main[model == m & term == t][[c]]
r <- function(s, c) rob[spec == s][[c]]
f1 <- data.table(
  lab = c("M2 Performance", "M3 Extended", "M5 Field (sample average)", "M7 Publisher (sample average)",
          "M4 Dynamic (short-run)", "CPD lagged three years"),
  est = c(g("M2", "lcpd", "est"), g("M3", "lcpd", "est"), num$avg5, num$avg7, g("M4", "lcpd", "est"),
          r("Lag-3 sample: CPD at t-3 (no window overlap)", "est")),
  se = c(g("M2", "lcpd", "se"), g("M3", "lcpd", "se"), num$avg5_se, num$avg7_se, g("M4", "lcpd", "se"),
         r("Lag-3 sample: CPD at t-3 (no window overlap)", "se")),
  kind = c(rep("Static specifications", 4), rep("Dynamic or longer lag", 2)))
f1[, lab := factor(lab, levels = rev(lab))]
p1 <- ggplot(f1, aes(est, lab, color = kind)) +
  geom_vline(xintercept = 0, color = INK, linewidth = 0.45) +
  geom_errorbarh(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), height = 0, linewidth = 0.5) +
  geom_point(size = 1.6) +
  geom_text(aes(label = sprintf("%.3f", est)), nudge_y = 0.33, size = 3.4, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("Static specifications" = BLUE, "Dynamic or longer lag" = RED),
                     breaks = c("Static specifications", "Dynamic or longer lag")) +
  scale_x_continuous(limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05)) +
  labs(x = "Elasticity of SJR with respect to CPD (95% CI)", y = NULL) + th +
  theme(axis.text = element_text(size = 10.5, color = INK), axis.title = element_text(size = 11, color = INK))
ggsave(file.path(FIG, "fig1_elasticity_across_models.png"), p1, width = W, height = 3.4, dpi = DPI, bg = "white")

## ---- Figure 2: field deviations (log CPD) ------------------------------
dev_fig <- function(tab, file, h) {
  t <- copy(tab); setorder(t, dev); t[, group := factor(group, levels = group)]
  t[, dir := factor(dir, levels = c("above", "below", "n.s."))]
  p <- ggplot(t, aes(dev, group, color = dir)) + geom_vline(xintercept = 0, color = INK, linewidth = 0.45) +
    geom_errorbarh(aes(xmin = dev - 1.96 * se, xmax = dev + 1.96 * se), height = 0, linewidth = 0.55) +
    geom_point(size = 1.7) + scale_color_manual(values = dircol, labels = dirlab, drop = FALSE) +
    labs(x = "Deviation from average elasticity (95% CI)", y = NULL) + th + inside_legend +
    theme(legend.background = element_blank(), plot.margin = margin(4, 8, 4, 4)) +
    guides(color = guide_legend(nrow = 2))
  ggsave(file.path(FIG, file), p, width = W, height = h, dpi = DPI, bg = "white")
}
## legend sits top-left in the field figure only if that corner is empty; put it bottom-right instead
f5 <- fread(O("lc_dev_m5.csv"))
dev_fig(f5, "fig2_field_deviations.png", 6.6)

## ---- Figure 3: the same data, opposite conclusions ---------------------
pp <- fread(O("lc_slope_paths_plus.csv"))
keep <- c("log(CPD), CPD > 0", "log(CPD + 0.1)", "asinh(CPD)", "log(1 + CPD)", "Implied elasticity from log(1 + CPD)")
pp <- pp[spec %in% keep]
pp[, idx := 100 * slope / slope[year == 2005], by = spec]
pp[, spec := factor(spec, levels = keep)]
cols <- c("log(CPD), CPD > 0" = BLUE, "log(CPD + 0.1)" = INK, "asinh(CPD)" = GREEN, "log(1 + CPD)" = RED,
          "Implied elasticity from log(1 + CPD)" = PURPLE)
ltys <- c("log(CPD), CPD > 0" = "solid", "log(CPD + 0.1)" = "solid", "asinh(CPD)" = "solid", "log(1 + CPD)" = "solid",
          "Implied elasticity from log(1 + CPD)" = "dashed")
pa <- ggplot(pp, aes(year, idx, color = spec, linetype = spec)) +
  geom_hline(yintercept = 100, color = INK, linewidth = 0.4) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.0) +
  scale_color_manual(values = cols) + scale_linetype_manual(values = ltys) +
  scale_x_continuous(breaks = c(2005, 2010, 2015, 2020, 2025)) +
  labs(x = NULL, y = "Slope on scaled CPD (2005 = 100)", title = "A. Slope on scaled CPD by year under four scalings") +
  guides(color = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2)) + th
## Panels B-C: same positive-CPD sample, observation-level group values, deviations from the sample average
gd <- fread(O("lc_group_devs.csv"))[level == "field" & sk %in% c("log", "l1p")]
w <- dcast(gd, group ~ sk, value.var = c("slope_dev", "el_dev"))   # slope_dev_l1p, el_dev_l1p, el_dev_log
w[, `:=`(dev_l1p = slope_dev_l1p, dev_log = el_dev_log, dev_c = el_dev_l1p)]
w[, short := group]; w[group == "Arts and Humanities", short := "Arts & Humanities"]
w[group == "Physics and Astronomy", short := "Physics"]
w[group == "Biochemistry, Genetics and Molecular Biology", short := "Biochemistry"]
w[group == "Immunology and Microbiology", short := "Immunology"]
rcor <- w[, cor(dev_log, dev_l1p)]
## labels are placed by hand (label x, label y, hjust) so that none overlaps a point or another label
sc_pts <- function(dt, xv, yv, lab_pos, ttl, xl, yl, diag = FALSE) {
  dt <- copy(dt); lp <- merge(lab_pos, dt[, .(group, short, px = get(xv), py = get(yv))], by = "group")
  p <- ggplot(dt, aes(.data[[xv]], .data[[yv]])) +
    geom_hline(yintercept = 0, color = INK, linewidth = 0.4) + geom_vline(xintercept = 0, color = INK, linewidth = 0.4)
  if (diag) p <- p + geom_abline(slope = 1, intercept = 0, color = RED, linewidth = 0.5, linetype = "dashed")
  else p <- p + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = RED, linewidth = 0.5, linetype = "dashed")
  p + geom_segment(data = lp, aes(x = px, y = py, xend = lx, yend = ly), inherit.aes = FALSE, color = INK, linewidth = 0.25) +
    geom_point(size = 1.8, color = BLUE) +
    geom_label(data = lp, aes(x = lx, y = ly, label = short, hjust = hj), inherit.aes = FALSE, size = 3.08, fontface = "bold",
               color = INK, fill = "white", label.size = 0, label.padding = unit(0.06, "lines")) +
    labs(x = xl, y = yl, title = ttl) + th + theme(plot.title = element_text(size = 11.5, face = "bold"))
}
LB <- data.table(group = c("Arts and Humanities", "Veterinary", "Chemistry", "Social Sciences", "Nursing", "Neuroscience",
                           "Biochemistry, Genetics and Molecular Biology"),
                 lx = c(0.10, -0.13, -0.235, 0.10, 0.205, 0.10, -0.15),
                 ly = c(-0.128, -0.088, -0.122, -0.068, -0.010, 0.110, 0.104),
                 hj = c(0.5, 0.5, 0, 0.5, 1, 0.5, 0.5))
pb <- sc_pts(w, "dev_l1p", "dev_log", LB, sprintf("B. Coefficients as estimated\n(r = %.2f)", rcor),
             "log(1 + CPD) coefficient, deviation", "log(CPD) elasticity, deviation")
cv <- w
rconv <- cv[, cor(dev_c, dev_log)]
LC <- data.table(group = LB$group,
                 lx = c(-0.17, -0.165, 0.02, -0.15, 0.07, 0.10, -0.03),
                 ly = c(-0.195, -0.075, -0.09, -0.015, -0.005, 0.135, 0.10),
                 hj = c(0, 0.5, 0, 0.5, 0, 0.5, 0.5))
pc <- sc_pts(cv, "dev_c", "dev_log", LC, sprintf("C. Converted to elasticities\n(r = %.2f)", rconv),
             "Converted elasticity, deviation", "log(CPD) elasticity, deviation", diag = TRUE)
lim <- c(-0.235, 0.165)
pc <- pc + coord_cartesian(xlim = lim, ylim = lim)
p3 <- pa / plot_spacer() / (pb + plot_spacer() + pc + plot_layout(widths = c(1, 0.06, 1))) + plot_layout(heights = c(1, 0.08, 1.15))
ggsave(file.path(FIG, "fig3_scaling.png"), p3, width = W, height = 7.8, dpi = DPI, bg = "white")
jsonlite::write_json(list(field_dev_cor_log_l1p = rcor, field_conv_cor = rconv), O("lc_fig_numbers.json"), auto_unbox = TRUE, digits = 10)

## ---- Figure 4: publisher deviations (log CPD) ---------------------------
d7 <- fread(O("lc_dev_m7.csv")); d7f <- fread(O("lc_dev_m7f.csv"))
t7 <- rbind(d7[, .(group, dev, se, dir, spec = "Publisher slopes only")], d7f[, .(group, dev, se, dir, spec = "Net of the 27 field slopes")])
t7[, group := factor(group, levels = d7[order(dev), group])]
t7[, spec := factor(spec, levels = c("Publisher slopes only", "Net of the 27 field slopes"))]
t7[, dir := factor(dir, levels = c("above", "below", "n.s."))]
p4 <- ggplot(t7, aes(dev, group, color = dir, shape = spec)) + geom_vline(xintercept = 0, color = INK, linewidth = 0.45) +
  geom_errorbarh(aes(xmin = dev - 1.96 * se, xmax = dev + 1.96 * se), height = 0, linewidth = 0.5, position = position_dodge(width = 0.72)) +
  geom_point(size = 1.7, position = position_dodge(width = 0.72)) +
  scale_color_manual(values = dircol, labels = dirlab, drop = FALSE) +
  scale_shape_manual(values = c("Publisher slopes only" = 16, "Net of the 27 field slopes" = 17)) +
  guides(color = guide_legend(order = 1, ncol = 1), shape = guide_legend(order = 2, ncol = 1)) +
  theme(legend.box = "horizontal") +
  labs(x = "Deviation from average elasticity (95% CI)", y = NULL) + th + inside_legend +
  theme(legend.background = element_blank(), legend.box = "horizontal", legend.box.just = "top")
ggsave(file.path(FIG, "fig4_publisher_deviations.png"), p4, width = W, height = 7.4, dpi = DPI, bg = "white")
cat("figures written to", FIG, "\n")
