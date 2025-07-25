

### plankton summary figures including vector plot for change in GCOB/GCOA and bubble plot for change in abundance values
### uses output/values from Plankton SDM file 


library(tidyr)

# Create df for GCOA data for zooplankton
distance_km_cala <- distance_meters_cala / 1000
distance_km_clad <- distance_meters_clad / 1000
distance_km_cyc  <- distance_meters_cyc / 1000
distance_km_larv <- distance_meters_larv / 1000

initial_bearing_cala <- bearing(center_hist_cala, center_proj_cala)
initial_bearing_clad <- bearing(center_hist_clad, center_proj_clad)
initial_bearing_cyc  <- bearing(center_hist_cyc, center_proj_cyc)
initial_bearing_larv <- bearing(center_hist_larv, center_proj_larv)

GCOB_Output <- data.frame(
  Taxa = c("Calanoids", "Cladocerans", "Cyclopoids", "Larvaceans"),
  Displacement_degrees = c(initial_bearing_cala, initial_bearing_clad, initial_bearing_cyc, initial_bearing_larv),
  Displacement_km = c(distance_km_cala, distance_km_clad, distance_km_cyc, distance_km_larv),
  stringsAsFactors = FALSE
)


# Now the fish larvae GCOA values
additional_GCOB <- data.frame(
  Taxa = c(
    "Billfish (Abiotic)",
    "Billfish (Biotic)",
    "Mahi-Mahi (Biotic)",
    "Mahi-Mahi (Abiotic)",
    "Other True Tuna (Biotic)",
    "Other True Tuna (Abiotic)",
    "Atl. Bluefin Tuna (Biotic)",
    "Atl. Bluefin Tuna (Abiotic)",
    "Skipjack Tuna (Biotic)",
    "Skipjack Tuna (Abiotic)",
    "Frigate Tuna (Biotic)",
    "Frigate Tuna (Abiotic)",
    "Little Tunny (Biotic)",
    "Little Tunny (Abiotic)"
  ),
  Displacement_degrees = c(
    initial_bearing_bill,
    initial_bearing_bill2,
    initial_bearing_cory2,
    initial_bearing_cory,
    initial_bearing_thunnus2,
    initial_bearing_thunnus,
    initial_bearing_thynnus2,
    initial_bearing_thynnus,
    initial_bearing_kat2,
    initial_bearing_kat,
    initial_bearing_aux2,
    initial_bearing_aux,
    initial_bearing_eut2,
    initial_bearing_eut
  ),
  Displacement_km = c(
    distance_km_bill,
    distance_km_bill2,
    distance_km_cory2,
    distance_km_cory,
    distance_km_thunnus2,
    distance_km_thunnus,
    distance_km_thynnus2,
    distance_km_thynnus,
    distance_km_kat2,
    distance_km_kat,
    distance_km_aux2,
    distance_km_aux,
    distance_km_eut2,
    distance_km_eut
  ),
  stringsAsFactors = FALSE
)

# Combine with existing GCOB_Output and transform
GCOB_Output <- rbind(GCOB_Output, additional_GCOB)
GCOB_Output$Displacement_degrees_trans <- GCOB_Output$Displacement_degrees - 180



GCOB_Output_filtered <- GCOB_Output %>%
  filter(Displacement_km >= 20)

vectorplot<-ggplot(GCOB_Output_filtered, aes(x = Displacement_degrees_trans, y = Displacement_km, color = Taxa)) +
  coord_radial(start = -0.9 * pi, end = pi) +
  geom_segment(aes(y = 0, xend = Displacement_degrees_trans, yend = Displacement_km), 
               arrow = arrow(length = unit(0.4, "cm"), type = 'closed'), size = .7) +
  scale_x_continuous(limits = c(-360, 0), breaks = c(-293, -189, -87, 16), 
                     labels = c("West", "North", "East", "South")) +  
  scale_y_continuous(limits = c(0, 225), breaks = c(0,50,100,150,200,250),
                     labels = c(" ","50 km","100 km","150 km","200 km","250 km")) +  
theme_bw() +
  theme(panel.border = element_blank(),  # Remove panel border
        panel.grid.major.x = element_blank(),  # Remove radial grid lines
        panel.grid.minor.x = element_blank(),  # Remove minor radial grid lines
        panel.grid.major.y = element_line(color = "darkgray"),  # Keep circular grid lines as dark gray
        panel.grid.minor.y = element_blank(),  # Optionally remove minor circular grid lines
        text = element_text(size = 8.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0.03, .99),
          legend.justification = c(.7, 1.02))+ # Adjust justification to top left corner
  scale_color_manual(
    values = c(
      "Cladocerans" = "palegreen3", 
      "Frigate Tuna (Abiotic)" = "cyan3", 
      "Billfish (Abiotic)" = "dodgerblue4", 
      "Atl. Bluefin Tuna (Biotic)" = "khaki3", 
      "Other True Tuna (Abiotic)" = "deepskyblue", 
      "Little Tunny (Abiotic)" = "lightsteelblue", 
      "Little Tunny (Biotic)" = "gold2", 
      "Skipjack Tuna (Abiotic)" = "royalblue",
      "Mahi-Mahi (Abiotic)" = "skyblue3"
    ))
vectorplot

ggsave("vectorplot.png",vectorplot, dpi = 300, bg = "white",
       width = 1961,
       height = 1400,
       units = "px") 

#percent change in abundnace for each taxa

# Reshape the data to long format
projections.df.subset <- projections.df %>%
  dplyr::select(cyclopoid_percentchange, calanoid_percentchange, larvacean_percentchange, 
                cladoceran_percentchange, billfish_percentchange, thunnus_percentchange, 
                thynnus_percentchange, cory_percentchange, aux_percentchange, eut_percentchange, kat_percentchange,
                billfish_percentchange2, thunnus_percentchange2, 
                thynnus_percentchange2, cory_percentchange2, aux_percentchange2, eut_percentchange2, kat_percentchange2)

projections_long <- projections.df.subset %>%
  pivot_longer(cols = everything(), names_to = "taxa", values_to = "percent_change") %>%
  mutate(taxa = case_when(
    taxa == "cyclopoid_percentchange" ~ "Cyclopoid",
    taxa == "calanoid_percentchange" ~ "Calanoid",
    taxa == "larvacean_percentchange" ~ "Larvacean",
    taxa == "cladoceran_percentchange" ~ "Cladoceran",
    taxa == "billfish_percentchange" ~ "Istiophoridae (Abiotic)",
    taxa == "thunnus_percentchange" ~ "Thunnus spp. (Abiotic)",
    taxa == "thynnus_percentchange" ~ "T. thynnus (Abiotic)",
    taxa == "eut_percentchange" ~ "E. alletteratus (Abiotic)",
    taxa == "aux_percentchange" ~ "Auxis (Abiotic)",
    taxa == "kat_percentchange" ~ "K. pelamis (Abiotic)",
    taxa == "cory_percentchange" ~ "Coryphaena (Abiotic)",  
    taxa == "billfish_percentchange2" ~ "Istiophoridae (Biotic)",
    taxa == "thunnus_percentchange2" ~ "Thunnus spp. (Biotic)",
    taxa == "thynnus_percentchange2" ~ "T. thynnus (Biotic)",
    taxa == "eut_percentchange2" ~ "E. alletteratus (Biotic)",
    taxa == "aux_percentchange2" ~ "Auxis (Biotic)",
    taxa == "kat_percentchange2" ~ "K. pelamis (Biotic)",
    taxa == "cory_percentchange2" ~ "Coryphaena (Biotic)"
  ),
  ModelType = case_when(
    taxa %in% c("Cyclopoid", "Calanoid", "Larvacean", "Cladoceran") ~ "Zooplankton",
    taxa %in% c("Auxis (Biotic)", "E. alletteratus (Biotic)", "K. pelamis (Biotic)", "Istiophoridae (Biotic)", "Thunnus spp. (Biotic)", "T. thynnus (Biotic)", "Coryphaena (Biotic)") ~ "PreyFish",
    taxa %in% c("Istiophoridae (Abiotic)", "Thunnus spp. (Abiotic)", "T. thynnus (Abiotic)", "Coryphaena (Abiotic)", "Auxis (Abiotic)", "E. alletteratus (Abiotic)", "K. pelamis (Abiotic)") ~ "ClimateFish"
  ))

projections_long <- projections_long %>%
  filter(!is.na(taxa))

mean_values <- projections_long %>%
  group_by(taxa, ModelType) %>%
  summarize(mean_percent_change = mean(percent_change))

# Manually adjust mean values to +100 if they exceed +100 (turn off for second plot, turn on for first plot)
#mean_values <- mean_values %>%
# mutate(mean_percent_change = if_else(mean_percent_change > 100, 100, mean_percent_change))

# Define the "Zooplankton" group taxa
zooplankton_taxa <- c("Calanoid", "Cladoceran", "Cyclopoid", "Larvacean")

# Create a custom order: "Zooplankton" group first, then other taxa in alphabetical order
ordered_taxa <- c(
  zooplankton_taxa,  # Zooplankton group
  sort(setdiff(unique(projections_long$taxa), zooplankton_taxa)) # Other taxa sorted alphabetically
)

# Apply this order to the `taxa` variable in both datasets
projections_long$taxa <- factor(projections_long$taxa, levels = ordered_taxa)
mean_values$taxa <- factor(mean_values$taxa, levels = ordered_taxa)

# Updated ggplot code
percentchangeplot <- ggplot(projections_long, aes(x = percent_change, y = taxa)) +
  geom_point(aes(size = 1, color = "gray99"), alpha = 0.3, show.legend = FALSE) + # Individual observations
  geom_point(data = mean_values, aes(x = mean_percent_change, y = taxa, color = ModelType),
             size = 7, alpha = 0.8, stroke = 2, show.legend = FALSE) + # Mean points with larger size
  scale_color_manual(values = c("Zooplankton" = "seagreen3", "PreyFish" = "khaki2", "ClimateFish" = "dodgerblue2")) + # Custom colors
  geom_vline(xintercept = 0, color = "black", size = 1) + # Vertical line at x = 0
  xlim(-100, 100) + # Set x-axis bounds
  theme_minimal() +
  labs(x = "Percent Change in Abundance", y = "Plankton Taxa") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 9)) +
  scale_y_discrete(labels = c(
    "Cladoceran" = expression(paste("Cladocerans")),
    "Calanoid" = expression(paste("Calanoids")),
    "Auxis (Abiotic)" = expression(paste("Frigate Tuna (Abiotic)")),
    "Auxis (Biotic)" = expression(paste("Frigate Tuna (Biotic)")),
    "Cyclopoid" = expression(paste("Cyclopoids")),
    "K. pelamis (Abiotic)" = expression(paste("Skipjack Tuna (Abiotic)")),
    "K. pelamis (Biotic)" = expression(paste("Skipjack Tuna (Biotic)")),
    "Larvacean" = expression(paste("Larvaceans")),
    "T. thynnus (Abiotic)" = expression(paste("Atl. Bluefin Tuna (Abiotic)")),
    "T. thynnus (Biotic)" = expression(paste("Atl. Bluefin Tuna (Biotic)")),
    "E. alletteratus (Abiotic)" = expression(paste("Little Tunny (Abiotic)")),
    "E. alletteratus (Biotic)" = expression(paste("Little Tunny (Biotic)")),
    "Coryphaena (Abiotic)" = expression(paste("Mahi-Mahi (Abiotic)")),
    "Coryphaena (Biotic)" = expression(paste("Mahi-Mahi (Biotic)")),
    "Istiophoridae (Abiotic)" = expression(paste("Billfish (Abiotic)")),
    "Istiophoridae (Biotic)" = expression(paste("Billfish (Biotic)")),
    "Thunnus spp. (Abiotic)" = expression(paste("Other True Tuna (Abiotic)")),
    "Thunnus spp. (Biotic)" = expression(paste("Other True Tuna (Biotic)"))
  ))

percentchangeplot

ggsave("percentchangeplot.png",percentchangeplot, dpi = 350, bg = "white",
       width = 1961,
       height = 2200,
       units = "px") 

### now make seperate plot for those taxa with really big values

# List of taxa with mean values > 100%
high_percent_taxa <- c("Coryphaena (Abiotic)", "Thunnus spp. (Abiotic)", 
                       "Istiophoridae (Abiotic)", "K. pelamis (Abiotic)")

# Custom y-axis labels for each taxon
custom_y_labels <- c(
  "Coryphaena (Abiotic)" = expression(paste("Mahi-Mahi (Abiotic)")),
  "Thunnus spp. (Abiotic)" = expression(paste("Other True Tuna (Abiotic)")),
  "Istiophoridae (Abiotic)" = expression(paste("Billfish (Abiotic)")),
  "K. pelamis (Abiotic)" = expression(paste("Skipjack Tuna (Abiotic)"))
)

# Create individual plots for each taxon
individual_plots <- lapply(1:length(high_percent_taxa), function(i) {
  taxon_name <- high_percent_taxa[i]
  
  # Filter data for the specific taxon
  taxon_data <- projections_long %>% filter(taxa == taxon_name)
  taxon_mean <- mean_values %>% filter(taxa == taxon_name)
  
  # Create the plot
  plot <- ggplot(taxon_data, aes(x = percent_change, y = taxa)) +
    geom_point(aes(size = 1, color = "gray99"), alpha = 0.3, show.legend = FALSE) + # Individual observations
    geom_point(data = taxon_mean, 
               aes(x = mean_percent_change, y = taxa, color = ModelType), 
               size = 7, alpha = 0.8, stroke = 2, show.legend = FALSE) + # Mean points
    scale_color_manual(values = c("Zooplankton" = "seagreen3", 
                                  "PreyFish" = "khaki2", 
                                  "ClimateFish" = "dodgerblue3")) + # Custom colors
    theme_minimal() +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 14))
  
  # Modify the y-axis labels for the specific taxon
  plot <- plot + scale_y_discrete(labels = custom_y_labels[taxon_name])
  
  # Remove the y-axis label
  plot <- plot + theme(axis.title.y = element_blank())  # Remove y-axis label
  
  # Add custom x-axis label only to the bottom plot
  if (i == length(high_percent_taxa)) {
    plot <- plot + labs(x = "Percent Change in Abundance")
  } else {
    plot <- plot + theme(axis.title.x = element_blank()) # Remove x-axis title for all but the bottom plot
  }
  
  # Set dynamic x-limits for each plot
  plot <- plot + xlim(range(taxon_data$percent_change, na.rm = TRUE))
  
  return(plot)
})

# Combine the plots using grid.arrange (stacking vertically)
highpercentchangeplot <- grid.arrange(grobs = individual_plots, ncol = 1)

ggsave("highpercentchangeplot.png",highpercentchangeplot, dpi = 250, bg = "white",
       width = 2200,
       height = 2000,
       units = "px") 



