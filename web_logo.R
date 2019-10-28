library(tidyverse)
library(rlang)
library(poissoned)
library(deldir)
library(rlist)
library(magrittr)
library(polyclip)
library(mgcv)

#takes a dataframe of polygons (with id column) and dataframe of points (x, y)
#and returns a dataframe of points that are inside the polygons and classified
#by the id of the original polygon. Format (x, y, id)
texturize <- function(polys, pts) {
  polys %>%
    select(x, y) %>%
    split(., polys$id) %>%
    map(., ~in.out(as.matrix(.x), as.matrix(pts))) %>%
    map(., ~cbind(.x, pts)) %>%
    map(., ~rename(.x, "inout" = ".x")) %>%
    map(., ~filter(.x, inout == TRUE)) %>%
    bind_rows(.id = "id") %>%
    select(-inout)
}

#functions for reformatting deldir outputs
cleanup <- function(x) x[ !names(x) %in% c("pt", "ptNum", "area", "id")]
cleanup2 <- function(x) x[ !names(x) %in% c("x", "y", "ptNum", "area", "id", "bp")]

#takes an input dataframe of seed points (x, y) and returns the voronoi tessellation as 
#a dataframe with (x, y, id) ready for ggplotting
voronize <- function(data) {
  tess <- deldir(as.data.frame(data))
  vor_list <- tile.list(tess) 
  
  vor_list %>%
    map( ~cleanup(.x)) %>%
    bind_rows(.id = "id") %>%
    select(id, x, y) %>%
    mutate(id = as.numeric(id))
}

#takes a set of points, and returns the points with their associated voronoi polygon IDs
#only reason this exists is to correlate polygon ID with seed point ID
get_seeds <- function(data) {
  tess <- deldir(as.data.frame(data))
  vor_list <- tile.list(tess) 
  
  vor_list %>%
    map( ~as.data.frame(t(unlist(cleanup2(.x))))) %>%
    bind_rows(.id = "id") %>%
    select(id, pt.x, pt.y) %>%
    mutate(id = as.numeric(id))
}

#takes a dataframe of polygons and applies a shape manipulation 
#polyclip::polyoffset has details on the parameters
shapify <- function(data, delta, jointype, miterlim = 2) {
  x_new <- split(data$x, data$id)
  y_new <- split(data$y, data$id)
  polygons <- Map(list, x = x_new, y = y_new)
  
  polygons2 <- lapply(polygons, polyoffset, delta = delta,
                      jointype = jointype, miterlim = miterlim)
  
  polygons2 %>%
    map(~as.data.frame(.x)) %>%
    bind_rows(.id = "id")
}


#draws num_points in a circle, badly
#min_r and max_r are the min and max radius possible
#jitter is the x y randomness
bad_circle <- function(num_pts = 10, min_r, max_r, jitter) {
  tibble(angle = seq(0, 2*pi, length.out = num_pts), r = sample(seq(min_r, max_r, length.out = 100), num_pts, replace = TRUE)) %>%
    mutate(x_jitter = sample(seq(-jitter, jitter, length.out = 100), num_pts, replace = TRUE), 
           y_jitter = sample(seq(-jitter, jitter, length.out = 100), num_pts, replace = TRUE),
           x = r*cos(angle) + x_jitter, 
           y = r*sin(angle) + y_jitter) %>%
    select(x, y) 
}

#draw concentric circles of pebbles then select some to keep
circle_pebbles <- function(num_pts = 10, min_r1 = 50, max_r1 = 100, jitter_1 = 100, min_r2 = 200, max_r2 = 300, jitter_2 = 100,
                           min_r3 = 350, max_r3 = 400, jitter_3 = 100, min_r4 = 500, max_r4 = 600, jitter_4 = 100,
                           min_r5 = 650, max_r5 = 700, jitter_5 = 100, expand = -30, round = 20, num_keepers = 8, probs = NULL) {
  circle1 <- bad_circle(num_pts, min_r1, max_r1, jitter_1)
  circle2 <- bad_circle(num_pts, min_r2, max_r2, jitter_2)
  circle3 <- bad_circle(num_pts, min_r3, max_r3, jitter_3)
  circle4 <- bad_circle(num_pts, min_r4, max_r4, jitter_4)
  circle5 <- bad_circle(num_pts, min_r5, max_r5, jitter_5)
  
  all_circles <- rbind(circle1, circle2, circle3, circle4, circle5)
  
  circular_layer <- voronize(all_circles) %>%
    shapify(delta = expand, jointype = "miter", miterlim = 1000) %>%
    shapify(delta = round, jointype = "round")
  
  keepers <- sample(1:50, num_keepers, prob = probs)
  
  seeds <- get_seeds(all_circles) %>%
    filter(id %in% keepers)
  
  pebbles <- 
    circular_layer %>%
    filter(id %in% keepers)
  
  list(seeds = seeds, pebbles = pebbles)
}



probs2 <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
# pebble_colors <- c('#E32C10', '#E35719', '#F9D3CE', '#133267')
# tex_colors <- c('#4d3d9a', '#f76975', '#ffffff', '#FCFBFA')

pebble_colors <- c('#20342a', '#f74713', '#686d2c', '#F9D3CE')
tex_colors <- c('#4d3d9a', '#E32C10', '#ffffff', '#eff0dd')
pts <- poisson_disc(ncols = 400, nrows = 400, cell_size = 5.5, k = 10, verbose = TRUE)
pts$x <- pts$x - 1000
pts$y <- pts$y - 1000

circ_layer5 <- circle_pebbles(probs = probs2, num_keepers = 20)[["pebbles"]] 

circ_layer6 <- 
  circ_layer5 %>%
  group_by(id) %>%
  group_map( ~mutate(.x, color = sample(pebble_colors, 1))) %>%
  ungroup()

tex_colored <- 
  texturize(circ_layer6, pts)  %>%
  mutate(size = sample(seq(0.03, 0.6, length.out = 100), size = nrow(.), replace = TRUE)) %>%
  group_by(id) %>%
  group_map( ~mutate(.x, color = sample(tex_colors, 1))) %>%
  ungroup()

tex_random <- 
  tex_colored %>%
  sample_n(nrow(tex_colored) / 10) %>%
  mutate(color = sample(tex_colors, nrow(.), replace = TRUE))

tex_final <- left_join(tex_colored, tex_random, by = c("x", "y")) %>%
  mutate(color.x = ifelse(is.na(color.y), color.x, color.y)) %>%
  select(id = id.x, x, y, size = size.x, color = color.x)

ggplot() +
  geom_polygon(data = circ_layer6, aes(x = x, y = y, group = id, fill = color)) +
  geom_point(data = tex_final, aes(x = x, y = y, size = size, color = color)) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_size_identity() +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", color = "NA"))

ggsave("header_pebbles2.png", height = 15, width = 15, bg = "transparent")



