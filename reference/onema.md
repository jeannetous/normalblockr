# French stream fish community data (ONEMA / OFB electrofishing surveys)

Fish biomass per species, aggregated by sampling station, together with
environmental covariates for each station. Derived from the fish
community monitoring of stream sections across metropolitan France
(1995-2018) by the French Office for Biodiversity (formerly ONEMA,
"Office National de l'Eau et des Milieux Aquatiques"), using
electrofishing.

## Usage

``` r
onema
```

## Format

A list with 2 elements:

- biomass:

  a 399 x 46 numeric matrix of total biomass (grams) per species
  (3-letter species code, columns) and station (row names).

- covariates:

  a data frame with one row per station, in the same order as the rows
  of \`biomass\`: \`station\` (matches \`rownames(biomass)\`) and 14
  environmental variables – \`slope\`, \`alt\` (altitude), \`d_source\`
  (distance to source), \`strahler\` (Strahler stream order),
  \`width_river_mean\`/\`width_river_cv\`, \`avg_depth_station_mean\`/
  \`avg_depth_station_cv\`, \`DBO_med\`/\`DBO_cv\` (biological oxygen
  demand), \`flow_med\`/\`flow_cv\` and
  \`temperature_med\`/\`temperature_cv\` (\`\_med\` = median, \`\_cv\` =
  coefficient of variation, over the monitoring period).

## Source

Danet, A., Mouchet, M., Bonnaffe, W., Thebault, E., Fontaine, C. (2021)
"Species richness and food-web structure jointly drive community biomass
and its temporal stability in fish communities", data set, Zenodo,
[doi:10.5281/zenodo.5095656](https://doi.org/10.5281/zenodo.5095656) .

## Details

The source data is organized around fishing \*operations\* (\`opcod\`):
several electrofishing operations can be carried out at the same
\*station\* over time. \`onema\` aggregates this to one row per station:
every operation is first mapped to its station
(\`fishing_protocol.csv\`), then biomass (in grams) is summed over every
operation recorded at a given station, separately for each species (i.e.
\`biomass\[s, \]\` is the total biomass of each species ever caught at
station \`s\`, not a single operation's catch). Stations without
environmental data are dropped. Species observed at fewer than 10
stations are also dropped: at that level of rarity, a zero-inflated
fit's iterative initialization can fit a species' handful of non-zero
observations exactly, driving its estimated noise precision to infinity
(these species also carry essentially no information for clustering
anyway).

## Examples

``` r
Y <- log(1 + onema$biomass)
X <- model.matrix(~ 1, data = onema$covariates)
nb_data <- NormalBlockData$new(Y, X)
# \donttest{
out <- normal_block(nb_data, 2:15, control = NB_control(clustering_init = "ward2"))
#> Fitting a diagonal normal-block-var model with unknown q 
#>   number of blocks = 2                number of blocks = 3                number of blocks = 4                number of blocks = 5                number of blocks = 6                number of blocks = 7                number of blocks = 8                number of blocks = 9                number of blocks = 10               number of blocks = 11               number of blocks = 12               number of blocks = 13               number of blocks = 14               number of blocks = 15           
#> DONE
# }
```
