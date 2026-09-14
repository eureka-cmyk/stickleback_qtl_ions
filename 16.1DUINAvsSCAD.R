library(tidyverse)
library(readxl)  

options(bitmapType = "cairo")

# --- 2. Data Import & Cleaning ---
file_path <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/DUINAvsSCAD.xlsx"

# Read the first row to get the correct element names
col_names <- names(read_excel(file_path, n_max = 0))

# Skip the first 2 rows (header + units) and apply the column names
# inomics_clean <- read_excel(file_path, skip = 2, col_names = col_names) %>%
#   # Remove potential empty rows at the bottom of Excel
#   filter(!is.na(ID))

inomics_clean <- read_excel(file_path, skip = 2, col_names = col_names) %>%
  filter(!is.na(ID)) %>%
  filter(POP %in% c("DUINA", "SCAD"))

# --- 3. Transformation for Distribution Plot ---
inomics_long <- inomics_clean %>%
  pivot_longer(cols = -(1:6), names_to = "element", values_to = "mg_kg") %>%
  # Filter out zeros or negative values for log transformation
  filter(mg_kg > 0)

# --- 4. Plotting: Elemental Distribution ---
# Using 'sans' family to avoid X11 font errors on the server
p1 <- ggplot(inomics_long, aes(x = element, y = mg_kg, color = POP)) +
  geom_jitter(alpha = 0.5, width = 0.2) +
  scale_y_log10() +
  theme_bw(base_family = "sans") + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  labs(title = "Elemental Distribution: DUINA vs SCAD", 
       y = "Concentration (mg/kg, log10 scale)")

p1



# Save the plot since X11 display might be disabled
ggsave("ionome_distribution_clean.png", p1, width = 12, height = 7)

# --- 5. Multi-Dimensional Scaling (MDS) ---
# Select only elemental columns (from 7th column onwards)
element_matrix <- inomics_clean %>% select(-(1:6))

# Standardize the data (Scale) is CRITICAL for ionomics
# This prevents major elements (like P) from dominating trace elements (like Cd)

dist_matrix <- dist(as.matrix(element_matrix))
dist_matrix <- dist(scale(element_matrix))
mds_res <- cmdscale(dist_matrix, k = 2)

# Combine MDS coordinates with metadata
inomics_mds <- inomics_clean %>%
  select(1:6) %>%
  mutate(mds1 = mds_res[,1], mds2 = mds_res[,2])

# --- 6. Plotting: MDS ---
p2 <- ggplot(inomics_mds, aes(x = mds1, y = mds2, color = POP, shape = SEX)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw(base_family = "sans") +
  labs(title = "MDS Plot of Stickleback Ionomics",
       # subtitle = "Based on standardized elemental profiles",
       x = "MDS Coordinate 1", 
       y = "MDS Coordinate 2")

p2
# Save the MDS plot
ggsave("ionome_mds_clean.png", p2, width = 8, height = 6)



# --- PCA ---
# Select only elemental columns
element_matrix <- inomics_clean %>% select(-(1:6))
element_matrix <- inomics_all %>%
  select(-ID, -POP, -SEX, -`Ti semi-quant`,-'Ti')
colSums(!is.finite(as.matrix(element_matrix)))
# Standardize before PCA
pca_res <- prcomp(element_matrix, center = TRUE, scale. = TRUE)

# Percentage of variance explained
var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

# Combine PCA coordinates with metadata
inomics_pca <- inomics_all %>%
  select(1:6) %>%
  mutate(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2]
  )

# --- Plotting: PCA ---
p_pca <- ggplot(inomics_pca, aes(x = PC1, y = PC2, color = POP, shape = SEX)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw(base_family = "sans") +
  labs(
    title = "PCA Plot of Stickleback Ionomics",
    x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
  )

p_pca

# Save PCA plot
ggsave("/gpfs01/home/shxyt1/data/sticklebacks/pheno/ionome_pca_clean.png", p_pca, width = 8, height = 6)



# --- 1. Load Libraries ---
library(tidyverse)
library(readxl)
library(scales) # For non-scientific notation

# --- 2. Load and Clean Original Excel Data ---
file_excel <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/DUINAvsSCAD.xlsx"
col_names_excel <- names(read_excel(file_excel, n_max = 0))

data_excel <- read_excel(file_excel, skip = 2, col_names = col_names_excel) %>%
  filter(!is.na(ID)) %>%
  filter(POP %in% c("DUINA", "SCAD")) %>%
  select(ID, POP, SEX, Na:U) # Select metadata and elements

# --- 3. Load and Clean New POND Data (CSV) ---
file_csv <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/qtl_all_phenotype_data2.csv"
data_pond <- read_csv(file_csv) %>%
  rename(ID = fish_id, SEX = sex) %>%
  mutate(POP = "POND") %>%
  # Select the same columns as Excel data to allow row binding
  select(ID, POP, SEX, Na:U)

# --- 4. Combine All Data ---
inomics_all <- bind_rows(data_excel, data_pond)

inomics_select <- inomics_all %>%
  select(ID, POP, SEX, Ca, Co, Fe, Sr, P, V)

# --- 5. Elemental Distribution Plot ---
inomics_long <- inomics_select %>%
  pivot_longer(cols = c(Ca, Co, Fe, Sr, P, V),
               names_to = "element",
               values_to = "mg_kg") %>%
  filter(mg_kg > 0)


p1 <- ggplot(inomics_long, aes(x = element, y = mg_kg, color = POP)) +
  # scale_y_log10() + 
  facet_wrap(~element, scale = "free") +
 # facet_wrap(~POP, nrow = 1) +
 # geom_jitter(position = position_dodge(width = 0.6), alpha = 0.4) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    alpha = 0.4
  )+
  theme_bw(base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 0), # Horizontal text
    legend.position = "right"
  ) +
  labs(title = "DUINA vs SCAD vs POND", y = "Concentration (mg/kg)")

p1
ggsave("/gpfs01/home/shxyt1/data/sticklebacks/pheno/DUINA vs SCAD vs PONDv3_indEle.png", p1, width = 10, height = 7)

p1.1 <- ggplot(inomics_long, aes(x = element, y = mg_kg, color = POP)) +
  scale_y_log10() + 
  # facet_wrap(~POP, nrow = 1) +
  # geom_jitter(position = position_dodge(width = 0.6), alpha = 0.4) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    alpha = 0.4
  )+
  theme_bw(base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 0), # Horizontal text
    legend.position = "right"
  ) +
  labs(title = "DUINA vs SCAD vs POND", y = "Concentration (mg/kg)")

p1.1
ggsave("/gpfs01/home/shxyt1/data/sticklebacks/pheno/DUINA vs SCAD vs PONDv3.png", p1, width = 10, height = 7)


# --- 1. Prepare Data ---
# Separate POND and others to control layers
data_pond_only <- inomics_long %>% filter(POP == "POND")
data_targets <- inomics_long %>% 
  filter(POP %in% c("DUINA", "SCAD")) %>%
  # Ensure the target points are plotted on top
  mutate(POP = factor(POP, levels = c("POND","DUINA","SCAD")))

# --- 2. Plotting ---
p_combined <- ggplot() +
  # Layer 1: POND as the Background Boxplot
  geom_boxplot(data = data_pond_only, 
               aes(x = element, y = mg_kg, fill = "POND"), 
               color = "green4", 
               outlier.shape = NA, width = 0.7) +
  
  # Layer 2: DUINA and SCAD as bright colored points
  # Using geom_jitter to spread points
  geom_jitter(data = data_targets, 
              aes(x = element, y = mg_kg, color = POP),
              width = 0.2,     
              height = 0,      
              # size = 2.5, 
              alpha = 0.5) +
  
  # 3. Adjust Scales & Labels
  scale_y_log10() +
  scale_color_discrete(name = "POP") +
  scale_fill_manual(name = NULL, values = c("POND" = "green")) +
  
  guides(
    color = guide_legend(order = 1), 
    fill = guide_legend(order = 2)  
  ) +
  
  theme_bw(base_family = "sans") +  
  # CRITICAL: Flip coordinates to make the plot vertical
  scale_fill_manual(name = NULL, values = c("POND" = "green")) +
  
  theme_bw(base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.y = element_text(size = 9), # Element names will be on the left
    # axis.title = element_text(face = "bold")
  ) +
  labs(title = "DUINA vs SCAD vs POND ",
       x = "Element", y = "Concentration (mg/kg, Log10)")

p_combined

ggsave("/gpfs01/home/shxyt1/data/sticklebacks/pheno/DUINA vs SCAD vs POND_V2.png", p_combined, width = 12, height = 7)




# --- 6. MDS Analysis with Scale ---
# Extract elements, ensure numeric, handle potential NAs
element_matrix <- inomics_all %>% 
  select(Na:U) %>% 
  mutate(across(everything(), as.numeric))

# Use Scale to give every element equal weight
dist_matrix <- dist(scale(element_matrix)) 
mds_res <- cmdscale(dist_matrix, k = 2)

inomics_mds <- inomics_all %>%
  select(ID, POP, SEX) %>%
  mutate(mds1 = mds_res[,1], mds2 = mds_res[,2])

# --- 7. MDS Plot ---
p2 <- ggplot(inomics_mds, aes(x = mds1, y = mds2, color = POP, shape = SEX)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("DUINA" = "#E41A1C", "SCAD" = "#377EB8", "POND" = "#4DAF4A")) +
  theme_bw(base_family = "sans") +
  labs(title = "MDS: Comparing Three Populations (Scaled)",
       x = "MDS1 (Standardized Unit)", y = "MDS2 (Standardized Unit)")

ggsave("mds_combined.png", p2, width = 9, height = 7)



#predict pond
train_data <- inomics_all %>% filter(POP %in% c("SCAD", "DUINA"))
test_data  <- inomics_all %>% filter(POP == "POND")

train_matrix <- train_data %>% select(Ca, Co, Fe, Sr, P, V)

pca_model <- prcomp(train_matrix, center = TRUE, scale. = TRUE)

test_matrix <- test_data %>% select(Ca, Co, Fe, Sr, P, V)

pond_proj <- predict(pca_model, newdata = test_matrix)

train_pca <- train_data %>%
  select(ID, POP, SEX) %>%
  mutate(
    PC1 = pca_model$x[,1],
    PC2 = pca_model$x[,2]
  )

pond_pca <- test_data %>%
  select(ID, POP, SEX) %>%
  mutate(
    PC1 = pond_proj[,1],
    PC2 = pond_proj[,2]
  )

pca_all <- bind_rows(train_pca, pond_pca)

p_predict <- ggplot(pca_all, aes(x = PC1, y = PC2, color = POP, shape = SEX)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw()


ggsave("/gpfs01/home/shxyt1/data/sticklebacks/pheno/predict.png", p_predict, width = 8, height = 6)
