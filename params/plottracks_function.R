# center a vector by dividing through a quantile (typically the median)
# use a pseudocount to prevent division by zero
.center_robust <- function(x, pscnt = 0.1, relquant = 0.5) {
    xc <- quantile(x + pscnt, prob = relquant)
    (x + pscnt) / xc
}

# smooth a vector using a running mean over windows of k_smooth elements
.smooth_runmean <- function(x, k_smooth) {
    require(S4Vectors, quietly = TRUE)
    as.numeric(runmean(x = as(x, "Rle"),
                       k = k_smooth, endrule = "constant"))
}

# plot genome-browser like tracks given a data.frame and gene names
plottracks <- function(df, gns, cols = NULL, k_smooth = 101,
                       lwd = 0.75, window_width = 13000, font_size = 14) {
    stopifnot(exprs = {
        require(tidyverse, quietly = TRUE)
        require(cowplot, quietly = TRUE)
        all(c("gene_name", "start", "end", "strand",
              "window_start", "window_end", "cpm_chip",
              "cpm_input") %in% colnames(df))
        window_width <= min(df |>
                                group_by(gene_name) |>
                                mutate(w = window_end - window_start + 1) |>
                                ungroup() |>
                                pull(w))
        all(gns %in% df$gene_name)
        is.null(cols) || all(gns %in% names(cols))
    })

    # define colours
    if (is.null(cols)) {
        cols <- hcl.colors(n = length(gns), palette = "Dark3")
        names(cols) <- gns
    }
    cols <- c(structure(unname(cols[gns]), names = paste0(gns, "_rel_chip")),
              structure(rep("gray", length(gns)), names = paste0(gns, "_rel_input")))


    # filter df, smooth and sub-sample
    pd <- df |>
        filter(gene_name %in% gns) |>
        group_by(gene_name) |>
        mutate(gene_name = factor(gene_name, levels = gns),
               pos = (seq(min(window_start), max(window_end), by = 1) -
                          ifelse(strand == "+", start, end)) *
                   ifelse(strand == "+", 1, -1),
               smooth_chip = .smooth_runmean(cpm_chip, k_smooth),
               smooth_input = .smooth_runmean(cpm_input, k_smooth),
               rel_chip = .center_robust(smooth_chip, pscnt = 0.001),
               rel_input = .center_robust(smooth_input, pscnt = 0.001)) |>
        filter(pos >= -(window_width / 2) & pos <= (window_width / 2)) |>
        slice(which(pos %% 40 == 0)) |>
        ungroup() |>
        pivot_longer(cols = starts_with("rel_"),
                     values_to = "rel",
                     names_to = "type") |>
        mutate(colour_type = paste0(gene_name, "_", type))

    # visualize
    pL <- lapply(split(seq.int(nrow(pd)), pd$gene_name), function(i) {
        nm <- as.character(pd$gene_name)[i[1]]
        ymax <- round(max(pd[i, "rel"]))
        xwidth <- pd[i[1], ] |> mutate(w = end - start + 1) |> pull(w)
        df0 <- data.frame(x = c(0, 0), y = c(0, 0.5) * ymax,
                          xend = c(0, xwidth), yend = c(0.5, 0.5) * ymax)
        p <- ggplot(pd[i, ], aes(pos, rel)) +
            geom_area(aes(fill = colour_type)) +
            geom_segment(data = df0[1, ], inherit.aes = FALSE, linewidth = lwd,
                         mapping = aes(x = x, xend = xend, y = y, yend = yend)) +
            geom_segment(data = df0[2, ], inherit.aes = FALSE, linewidth = lwd,
                         arrow = arrow(length = unit(0.075, "npc"),
                                       angle = 25, type = "closed"),
                         mapping = aes(x = x, xend = xend, y = y, yend = yend)) +
            geom_text(data = df0[2, ], inherit.aes = FALSE,
                      mapping = aes(x = x + 0.04 * window_width, y = y),
                      label = ifelse(grepl("^SP", nm),
                                     nm,
                                     paste0(nm, "\U207A")), # superscript '+'
                      fontface = ifelse(grepl("^SP", nm),
                                        "plain",
                                        "italic"),
                      hjust = 0, vjust = -0.4,
                      size = 3.88 / 14 * font_size * 0.9) +
            scale_fill_manual(values = cols) +
            scale_x_continuous(expand = expansion(mult = 0.01)) +
            scale_y_continuous(position = "right", breaks = c(0, ymax),
                               expand = expansion(mult = c(0, 0.01))) +
            coord_cartesian(clip = "off") +
            theme_cowplot(font_size = font_size) +
            theme(legend.position = "none",
                  axis.title.x = element_blank(),
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.y = element_text(size = 11 / 14 * font_size * 0.9)) +
            annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1,
                     label = .capitalize(nm), size = 3.88 / 14 * font_size)
        if (nm == gns[[1]]) {
            df0 <- data.frame(
                x = max(pd[i, "pos"]) - 2000, y = max(pd[i, "rel"]) * 0.75,
                xmid = max(pd[i, "pos"]) - 1000,
                xend = max(pd[i, "pos"]), yend = max(pd[i, "rel"]) * 0.75
            )
            p <- p +
                geom_segment(data = df0, inherit.aes = FALSE, linewidth = lwd,
                             mapping = aes(x = x, xend = xend, y = y, yend = yend)) +
                geom_text(data = df0, inherit.aes = FALSE, label = "2 kb",
                          mapping = aes(x = xmid, y = y),
                          hjust = 0.5, vjust = -0.4, size = 3.88 / 14 * font_size)
        }
        if (nm == "zas1") {
            p <- p + geom_text(
                data = data.frame(pos = -300, rel = max(pd[i, "rel"])),
                label = "*", size = 6 / 14 * font_size,
                hjust = 0.5, vjust = -0.2)
        }
        p
    })
    plot_grid(plotlist = pL, align = "v", ncol = 1)
}
