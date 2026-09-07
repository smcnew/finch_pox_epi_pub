# Plot results of PSEMS

library(dplyr)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)


# Load data

main_sem <- read.delim("output_tables/sem_environment_elevation.csv", sep = ",") %>%
    filter(P.Value < 0.05)  %>% # keep only significant paths
   mutate(stars = case_when(P.Value < 0.001 ~ "***",
                        P.Value < 0.01 ~ "**",
                        P.Value < 0.05 ~  "*",
                        TRUE ~ "")) %>%
  mutate(across(.cols = c(Response, Predictor), .fns =
           ~recode(.x,
                  prop_inf = "Pox prevalence",
                  scale_dens = "Host abundance",
                  scale_shannons = "Host diversity",
                  scale_templ2 = "Temperature",
                  scale_elev = "Elevation",
                  scale_precipl2 = "Precipitation"))) %>%
  mutate(style = ifelse(Std.Estimate < 0, "dashed", "solid"))

# PSEM results testing for effects of species abundance on pox prevalence

species <- read.delim("output_tables/species_psem.csv", sep = ",") %>%
  mutate(stars = case_when(P.Value < 0.001 ~ "***",
                           P.Value < 0.01 ~ "**",
                           P.Value < 0.05 ~  "*",
                           TRUE ~ "")) %>%
  mutate(across(.cols = c(Response, Predictor), .fns =
                  ~recode(.x,
                          prop_inf = "Pox prevalence",
                          scale_templ2 = "Temperature",
                          scale_elevation = "Elevation",
                          scale_precipl2 = "Precipitation",
                          scale_humidityl2 = "Humidity"))) %>%
  mutate(style = ifelse(P.Value < 0.05, "solid", "dotted"))

# Fuliginosa
ful <- species %>%
  filter(species == "FUL") %>%
  mutate(across(.cols = c(Response, Predictor), .fns =
                  ~recode(.x,
                          scale_prop_FUL = "Proportion G. fuliginosa")))
# Scandens
sca <- species %>%
  filter(species == "SCA") %>%
  mutate(across(.cols = c(Response, Predictor), .fns =
                  ~recode(.x,
                          scale_prop_SCA = "Proportion G. scandens")))
# Magnirostris
mag <- species %>%
  filter(species == "MAG") %>%
  mutate(across(.cols = c(Response, Predictor), .fns =
                  ~recode(.x,
                          scale_prop_MAG = "Proportion G. magnirostris")))
#
# Diagram main SEM --------------------------------------------------------


# Format labels
main_sem <- main_sem %>%
  mutate(Label = paste0(round(Std.Estimate, 3),
                       " [",
                       round(Lower_CI, 3),
                       " — ",
                       round(Upper_CI, 3),
                       "] ",
                       stars)) # Create a label

# Format edges

main_sem <- main_sem %>%
  mutate(edges = paste0(
    "  \"", main_sem$Predictor, "\" -> \"", main_sem$Response,
    "\" [label = \"", main_sem$Label,
    "\", style = ", main_sem$style,
    "]"
  ))


# Create node list
nodes <- unique(c(main_sem$Predictor, main_sem$Response))

nodes <- paste0(
  "  \"", nodes, "\" [shape = box]"
)

# Assemble Graphviz code
dot_code <- paste0(
  "digraph SEM {\n",
  "  graph [layout = dot, rankdir = TB]\n",
  "  node [fontname = Helvetica, fontsize = 12]\n",
  "  edge [fontname = Helvetica, fontsize = 11]\n\n",
 "{ rank = same; Temperature; Precipitation, Elevation}",
  paste(nodes, collapse = "\n"),
  "\n\n",
  paste(main_sem$edges, collapse = "\n"),
  "\n}"
)

# Create plot
path_diagram_main <- grViz(dot_code)

# Convert to SVG for export
svg <- export_svg(path_diagram_main)

# Write file
rsvg_pdf(charToRaw(svg),
             file = "output_plots/path_diagram_main_sem_revision.pdf")

#
# Path diagrams for G. fuliginosa -----------------------------------

# Format labels
ful <- ful %>%
  mutate(Label = paste0(round(Std.Estimate, 3),
                        " [",
                        round(Lower_CI, 3),
                        " — ",
                        round(Upper_CI, 3),
                        "] ",
                        stars)) # Create a label

# Format edges

ful <- ful %>%
  mutate(edges = paste0(
    "  \"", Predictor, "\" -> \"", Response,
    "\" [label = \"", Label,
    "\", style = ", style,
    "]"
  ))


# Create node list
nodes_ful <- unique(c(ful$Predictor, ful$Response))

nodes_ful <- paste0(
  "  \"", nodes_ful, "\" [shape = box]"
)

# Assemble Graphviz code
dot_code_ful <- paste0(
  "digraph SEM {\n",
  "  graph [layout = dot, rankdir = TB]\n",
  "  node [fontname = Helvetica, fontsize = 12]\n",
  "  edge [fontname = Helvetica, fontsize = 11]\n\n",
  "{ rank = same; Temperature; Precipitation, Elevation, Humidity}",
  paste(nodes_ful, collapse = "\n"),
  "\n\n",
  paste(ful$edges, collapse = "\n"),
  "\n}"
)

# Create plot
ful_plot <- grViz(dot_code_ful)


# Convert to SVG for export
svg_ful <- export_svg(ful_plot)

# Write file
rsvg_pdf(charToRaw(svg_ful),
         file = "output_plots/path_diagram_fuliginosa.pdf")



# Plot for G. scandens ----------------------------------------------------

# Format labels
sca <- sca %>%
  mutate(Label = paste0(round(Std.Estimate, 3),
                        " [",
                        round(Lower_CI, 3),
                        " — ",
                        round(Upper_CI, 3),
                        "] ",
                        stars)) # Create a label

# Format edges

sca <- sca %>%
  mutate(edges = paste0(
    "  \"", Predictor, "\" -> \"", Response,
    "\" [label = \"", Label,
    "\", style = ", style,
    "]"
  ))


# Create node list
nodes_sca <- unique(c(sca$Predictor, sca$Response))

nodes_sca <- paste0(
  "  \"", nodes_sca, "\" [shape = box]"
)

# Assemble Graphviz code
dot_code_sca <- paste0(
  "digraph SEM {\n",
  "  graph [layout = dot, rankdir = TB]\n",
  "  node [fontname = Helvetica, fontsize = 12]\n",
  "  edge [fontname = Helvetica, fontsize = 11]\n\n",
  "{ rank = same; Temperature; Precipitation, Elevation, Humidity}",
  paste(nodes_sca, collapse = "\n"),
  "\n\n",
  paste(sca$edges, collapse = "\n"),
  "\n}"
)

# Create plot
sca_plot <- grViz(dot_code_sca)


# Convert to SVG for export
svg_sca <- export_svg(sca_plot)

# Write file
rsvg_pdf(charToRaw(svg_sca),
         file = "output_plots/path_diagram_scandens.pdf")


# G. magnirostris ---------------------------------------------------------

# Format labels
mag <- mag %>%
  mutate(Label = paste0(round(Std.Estimate, 3),
                        " [",
                        round(Lower_CI, 3),
                        " — ",
                        round(Upper_CI, 3),
                        "] ",
                        stars)) # Create a label

# Format edges

mag <- mag %>%
  mutate(edges = paste0(
    "  \"", Predictor, "\" -> \"", Response,
    "\" [label = \"", Label,
    "\", style = ", style,
    "]"
  ))


# Create node list
nodes_mag <- unique(c(mag$Predictor, mag$Response))

nodes_mag <- paste0(
  "  \"", nodes_mag, "\" [shape = box]"
)

# Assemble Graphviz code
dot_code_mag <- paste0(
  "digraph SEM {\n",
  "  graph [layout = dot, rankdir = TB]\n",
  "  node [fontname = Helvetica, fontsize = 12]\n",
  "  edge [fontname = Helvetica, fontsize = 11]\n\n",
  "{ rank = same; Temperature; Precipitation, Elevation, Humidity}",
  paste(nodes_mag, collapse = "\n"),
  "\n\n",
  paste(mag$edges, collapse = "\n"),
  "\n}"
)

# Create plot
mag_plot <- grViz(dot_code_mag)


# Convert to SVG for export
svg_mag <- export_svg(mag_plot)

# Write file
rsvg_pdf(charToRaw(svg_mag),
         file = "output_plots/path_diagram_magnirostris.pdf")
