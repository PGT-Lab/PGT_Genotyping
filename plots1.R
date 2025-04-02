pacman::p_load(ggplot2, ggthemr, envalysis)

dados <- data.frame(
  Dado = c(
    "PC_Kings_2016_2017",
    "PI_ESALQ_2020",
    "PSC_Kings_2019",
    "PSC_Kings_2016_2019",
    "TP_USP_2020",
    "SC_CHOP_2017",
    "PS_CHOP_2017"),
  chip = c("Psych", "GSA", "GSA", "Psych", "GSA", "Omni", "Omni"),
  versions = c("v1.1", "v2.0", "v3.0", "v1.1", "v1.0", "v1.0", "v1.1"),
  Variantes = c(587111, 665469, 654014, 587111, 618495, 729858, 719024))
dados

ggthemr("light")

ggplot(dados, aes(Dado, Variantes, fill = chip)) +
geom_col() +
geom_text(label = dados$Variantes, hjust = 0.5, vjust = 1.5, color = "black", size = 2) +
geom_text(label = dados$versions, hjust = 0.5, vjust = 3, color = "black", size = 2) +
theme_publish(base_size = 7) +
theme(
 axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.7),
 axis.title.x = element_blank(),
 legend.position = "right")

ggsave("chip_variants.png", device = "png", height = 70, width = 100, units = "mm")

intersections <- data.frame(
  Dataset1 = c(
    "PC_Kings_2016_2017", "PC_Kings_2016_2017", "PC_Kings_2016_2017", 
    "PC_Kings_2016_2017", "PC_Kings_2016_2017", "PC_Kings_2016_2017", 
    "PI_ESALQ_2020", "PI_ESALQ_2020", "PI_ESALQ_2020", "PI_ESALQ_2020", 
    "PI_ESALQ_2020", "PSC_Kings_2019", "PSC_Kings_2019", "PSC_Kings_2019", 
    "PSC_Kings_2019", "PSC_Kings_2016_2019", "PSC_Kings_2016_2019", 
    "PSC_Kings_2016_2019", "TP_USP_2020", "TP_USP_2020", "SC_CHOP_2017"),
  Chip1 = c("Psych", "Psych", "Psych", "Psych", "Psych", "Psych", 
            "GSA", "GSA", "GSA", "GSA", "GSA", "GSA", "GSA", "GSA", "GSA", 
            "Psych", "Psych", "Psych", "GSA", "GSA", "Omni"),
  v1 = c("v1.1", "v1.1", "v1.1", "v1.1", "v1.1", "v1.1", 
               "v2.0", "v2.0", "v2.0", "v2.0", "v2.0", "v3.0", "v3.0", "v3.0", "v3.0", 
               "v1.1", "v1.1", "v1.1", "v1.0", "v1.0", "v1.0"),
  Dataset2 = c(
    "PI_ESALQ_2020", "PSC_Kings_2019", "PSC_Kings_2016_2019", "TP_USP_2020", 
    "SC_CHOP_2017", "PS_CHOP_2017", "PSC_Kings_2019", "PSC_Kings_2016_2019", 
    "TP_USP_2020", "SC_CHOP_2017", "PS_CHOP_2017", "PSC_Kings_2016_2019", 
    "TP_USP_2020", "SC_CHOP_2017", "PS_CHOP_2017", "TP_USP_2020", "SC_CHOP_2017", 
    "PS_CHOP_2017", "SC_CHOP_2017", "PS_CHOP_2017", "PS_CHOP_2017"),
  Chip2 = c("GSA", "GSA", "Psych", "GSA", "Omni", "Omni", "GSA", "Psych", 
            "GSA", "Omni", "Omni", "Psych", "GSA", "Omni", "Omni", "GSA", "Omni", 
            "Omni", "Omni", "Omni", "Omni"),
  v2 = c("v2.0", "v3.0", "v1.1", "v1.0", "v1.0", "v1.1", "v3.0", "v1.1", 
               "v1.0", "v1.0", "v1.1", "v1.1", "v1.0", "v1.0", "v1.1", "v1.0", 
               "v1.0", "v1.1", "v1.0", "v1.1", "v1.1"),
  Variantes = c(109542, 108461, 587111, 107784, 258365, 256052, 
                   635892, 109542, 613312, 139140, 138358, 108461, 589402, 137739, 
                   136991, 107784, 258365, 256052, 133277, 132530, 719023))

ggplot(intersections, aes(Dado, Variantes, fill = chip)) +
geom_col() +
geom_text(label = intersections$Variantes, hjust = 0.5, vjust = 1.5, color = "black", size = 2) +
geom_text(label = intersections$versions, hjust = 0.5, vjust = 3, color = "black", size = 2) +
theme_publish(base_size = 7) +
theme(
 axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.7),
 axis.title.x = element_blank(),
 legend.position = "right")

library(forcats)

# Prepare data
bar_df <- data.frame(
  intersection = names(intersection_data),
  value = intersection_data
)

pacman::p_load(ggalluvial)
library(ggplot2)
library(dplyr)
library(tidyr)
# Prepare the data from your table
intersection_df <- data.frame(
  Conjunto = c(
    "PC_Kings_2016_2017_Psychv1.1_Hg19 & PI_ESALQ_2020_GSAv2.0_Hg19",
    "PC_Kings_2016_2017_Psychv1.1_Hg19 & PSC_Kings_2019_GSAv3.0_Hg19",
    "PC_Kings_2016_2017_Psychv1.1_Hg19 & PSC_Kings_2016_2019_Psychv1.1_Hg19",
    "PC_Kings_2016_2017_Psychv1.1_Hg19 & TP_USP_2020_GSAv1.0_Hg19",
    "PC_Kings_2016_2017_Psychv1.1_Hg19 & SC_CHOP_2017_OmniExpressv1.0_Hg19",
    "PC_Kings_2016_2017_Psychv1.1_Hg19 & PS_CHOP_2017_OmniExpressv1.1_Hg19",
    "PI_ESALQ_2020_GSAv2.0_Hg19 & PSC_Kings_2019_GSAv3.0_Hg19",
    "PI_ESALQ_2020_GSAv2.0_Hg19 & PSC_Kings_2016_2019_Psychv1.1_Hg19",
    "PI_ESALQ_2020_GSAv2.0_Hg19 & TP_USP_2020_GSAv1.0_Hg19",
    "PI_ESALQ_2020_GSAv2.0_Hg19 & SC_CHOP_2017_OmniExpressv1.0_Hg19",
    "PI_ESALQ_2020_GSAv2.0_Hg19 & PS_CHOP_2017_OmniExpressv1.1_Hg19",
    "PSC_Kings_2019_GSAv3.0_Hg19 & PSC_Kings_2016_2019_Psychv1.1_Hg19",
    "PSC_Kings_2019_GSAv3.0_Hg19 & TP_USP_2020_GSAv1.0_Hg19",
    "PSC_Kings_2019_GSAv3.0_Hg19 & SC_CHOP_2017_OmniExpressv1.0_Hg19",
    "PSC_Kings_2019_GSAv3.0_Hg19 & PS_CHOP_2017_OmniExpressv1.1_Hg19",
    "PSC_Kings_2016_2019_Psychv1.1_Hg19 & TP_USP_2020_GSAv1.0_Hg19",
    "PSC_Kings_2016_2019_Psychv1.1_Hg19 & SC_CHOP_2017_OmniExpressv1.0_Hg19",
    "PSC_Kings_2016_2019_Psychv1.1_Hg19 & PS_CHOP_2017_OmniExpressv1.1_Hg19",
    "TP_USP_2020_GSAv1.0_Hg19 & SC_CHOP_2017_OmniExpressv1.0_Hg19",
    "TP_USP_2020_GSAv1.0_Hg19 & PS_CHOP_2017_OmniExpressv1.1_Hg19",
    "SC_CHOP_2017_OmniExpressv1.0_Hg19 & PS_CHOP_2017_OmniExpressv1.1_Hg19"
  ),
  Interseção = c(
    109542, 108461, 587111, 107784, 258365, 256052,
    635892, 109542, 613312, 139140, 138358,
    108461, 589402, 137739, 136991,
    107784, 258365, 256052,
    133277, 132530, 719023
  )
)

# Split the Conjunto column into two separate columns
plot_data <- intersection_df %>%
  separate(Conjunto, into = c("From", "To"), sep = " & ") %>%
  mutate(From = trimws(From), To = trimws(To))

# Create the alluvial plot
ggplot(plot_data,
       aes(axis1 = From, axis2 = To, y = Interseção)) +
  geom_alluvium(aes(fill = From), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/6, fill = "gray90", color = "gray50") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3,
            min.y = 10000) + # Only show labels for flows > 10,000
  scale_x_discrete(limits = c("From", "To"),
                   expand = c(0.05, 0.05)) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Dataset Intersections",
       y = "Intersection Size",
       x = "Dataset Relationships") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        panel.grid = element_blank())
