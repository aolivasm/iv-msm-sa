labels <- c("italic(a)==0.55*\",\"~~italic(b)==0.075", "italic(a)==1.10*\",\"~~italic(b)==0.075", "italic(a)==1.65*\",\"~~italic(b)==0.075", 
    "italic(a)==1.10*\",\"~~italic(b)==0.150")

decorate <- function(d) {
    d$outcome <- factor(d$outcome, levels = c("Continuous", "Binary"))
    d$setting_label <- factor(labels[d$setting], levels = labels)
    if ("display" %in% names(d)) 
        d$display <- factor(d$display, levels = c("IV-MSM", "Tan", "Dorn"))
    d
}

main <- decorate(main)
ep <- decorate(ep)
palette <- old$application_palette

theme_paper <- theme_bw(base_size = 10, base_family = "sans") + theme(panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.25), strip.background = element_rect(fill = "grey94", 
        colour = "grey70", linewidth = 0.3), strip.text = element_text(size = 10), legend.position = "bottom", 
    legend.title = element_blank(), plot.title = element_text(size = 12, face = "bold"), plot.caption = element_text(size = 9, 
        hjust = 0, lineheight = 1.1), panel.spacing = grid::unit(0.55, "lines"), plot.margin = margin(8, 
        10, 8, 8))

facetlab <- labeller(setting_label = label_parsed, outcome = label_value)

facet <- facet_grid(outcome ~ setting_label, labeller = facetlab)

method_scales <- list(scale_colour_manual(values = palette, drop = FALSE), scale_shape_manual(values = c(`IV-MSM` = 16, 
    Tan = 17, Dorn = 15)), scale_linetype_manual(values = c(`IV-MSM` = "solid", Tan = "dashed", Dorn = "dotted")))

ann <- data.frame(outcome = factor("Continuous", levels = c("Continuous", "Binary")), setting_label = factor(labels[4], 
    levels = labels), Lambda = 1.76)

pA <- ggplot(main, aes(Lambda, 100 * outward, colour = display, shape = display, linetype = display, 
    group = display)) + geom_hline(yintercept = 95, colour = "#B2182B", linetype = 2, linewidth = 0.4) + 
    geom_line(linewidth = 0.45, position = position_dodge(0.05)) + geom_point(size = 1.8, position = position_dodge(0.05)) + 
    geom_errorbar(aes(ymin = 100 * mc_lower, ymax = 100 * mc_upper), width = 0.025, linewidth = 0.35, 
        position = position_dodge(0.05)) + facet + method_scales + scale_x_continuous(breaks = c(1.8, 
    2, 2.5, 3), limits = c(1.72, 3.08)) + scale_y_continuous(limits = c(85, 101), breaks = c(85, 90, 
    95, 100)) + geom_text(data = ann, aes(x = Lambda, y = 86.3, label = "IV-MSM omitted\nfor incompatible values"), 
    inherit.aes = FALSE, hjust = 0, size = 2.55) + labs(title = "A. Outward confidence-set coverage", 
    x = NULL, y = "Coverage (%)") + theme_paper + guides(colour = "none", shape = "none", linetype = "none")

width_limits <- decorate(data.frame(outcome = rep(c("Continuous", "Binary"), each = 2), setting = 1, 
    Lambda = 2, confidence_width = c(0, 1.8, 0, 0.72)))

pB <- ggplot(main, aes(Lambda, confidence_width, colour = display, shape = display, linetype = display, 
    group = display)) + geom_blank(data = width_limits, aes(Lambda, confidence_width), inherit.aes = FALSE) + 
    geom_line(linewidth = 0.5) + geom_point(size = 1.8) + facet_grid(outcome ~ setting_label, scales = "free_y", 
    labeller = facetlab) + method_scales + scale_x_continuous(breaks = c(1.8, 2, 2.5, 3), limits = c(1.72, 
    3.08)) + scale_y_continuous(expand = expansion(mult = c(0, 0.025))) + geom_text(data = ann, aes(x = Lambda, 
    y = 0.22, label = "IV-MSM omitted\nfor incompatible values"), inherit.aes = FALSE, hjust = 0, size = 2.55) + 
    labs(title = "B. Mean outward confidence-set width", x = expression(Lambda), y = "Mean width") + 
    theme_paper

long <- rbind(transform(ep, endpoint = "Lower", estimate = lower, target = target_lower), transform(ep, 
    endpoint = "Upper", estimate = upper, target = target_upper))

long$x <- match(long$Lambda, c(1.8, 2, 2.5, 3)) + c(`IV-MSM` = -0.23, Tan = 0, Dorn = 0.23)[as.character(long$display)]

refs <- unique(long[c("outcome", "setting_label", "Lambda", "display", "endpoint", "x", "target")])

refs <- refs[is.finite(refs$target), ]

counts <- table(interaction(long$outcome, long$setting, long$Lambda, long$display, long$endpoint, drop = TRUE))

fig2 <- ggplot(long, aes(x, estimate, group = interaction(Lambda, display, endpoint), fill = display, 
    colour = display)) + geom_boxplot(position = "identity", width = 0.18, linewidth = 0.28, outlier.size = 0.35, 
    outlier.alpha = 0.45, alpha = 0.18) + geom_segment(data = refs, aes(x = x - 0.105, xend = x + 0.105, 
    y = target, yend = target), inherit.aes = FALSE, colour = "#C51B2D", linewidth = 0.6) + facet_grid(outcome ~ 
    setting_label, scales = "free_y", labeller = facetlab) + scale_colour_manual(values = palette) + 
    scale_fill_manual(values = palette) + scale_x_continuous(breaks = 1:4, labels = c("1.8", "2", "2.5", 
    "3"), limits = c(0.5, 4.5)) + labs(x = expression(Lambda), y = "ATE endpoint estimate", caption = paste("Stress designs; n = 6,000; 500 repetitions per box. Lower and upper endpoints share each method column.", 
    "Red segments: method-specific population endpoints. All outliers are shown; scales are shared within outcome rows.", 
    "Setting 4 includes incompatible IV-MSM estimates at Lambda = 2, without population references. Tan and Dorn adjust for X.", 
    sep = "\n")) + theme_paper + guides(colour = "none", fill = guide_legend(override.aes = list(alpha = 0.5)))

