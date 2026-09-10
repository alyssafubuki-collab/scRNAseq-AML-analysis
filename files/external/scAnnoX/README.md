# scAnnoX - An R Package Integrating Multiple Public Tools for Single-Cell Annotation
scAnnoX is an R package that integrates 10 different cell identity algorithms for single-cell sequencing data into a unified framework, and allows comparison between different algorithms.

# Install
To install scAnnoX, we recommend using "devtools":
```R
#install.packages("devtools")  
devtools::install_github('XQ-hub/scAnnoX')
```

# Prepare input data
```R
library(scAnnoX)
library(SingleCellExperiment)

# Import single cell profiles.
test.obj <- system.file('data', 'test.obj.rds', package = 'scAnnoX') %>% readRDS(.)
ref.obj <- system.file('data', 'ref.obj.rds', package = 'scAnnoX') %>% readRDS(.)
```

# Set parameters
| Parameters   | Description |
| ------------ | ------------------------------------------ |
| obj.seu      | Seurat object, which needs to be annotated. |
| ref.obj      | Seurat object, only when used with reference-based tools. Default: NULL.  |
| ref.ctype    | The cell type column in ref.obj, only when used with reference-based tools. Default: NULL.   |
| marker.lst   | A list contained maker genes for each cell type for marker based tool.   |
| method       | A vector of automated annotation tools.   |
| select.marker| Specify method for inferring markers for each subset. Default: Seurat.  |
| top.k        | Top k expressed genes of each subset remained. Default: NULL.  |
| strategy     | Category of single cell annotation tool.   |

# Single cell annotation via scAnnoX
Exclusively employed a singular single-cell annotation tool.
```R
# Take one of the reference-based annotation tools for example.
pred.obj <- autoAnnoTools(
    test.obj,
    ref.obj = ref.obj,
    ref.ctype = 'CellType',
    marker.lst = marker.lst,
    method = 'SingleR',
    select.marker = 'Seurat',
    top.k = 30,
    strategy = 'reference-based'
) 

# Take one of the marker-based annotation tools for example.
Idents(ref.obj) <- ref.obj$CellType
marker.lst <- findMarkerToolsForSc(ref.obj, to.list = TRUE, top.k = 30)
pred.obj <- autoAnnoTools(
    test.obj,
    ref.obj = ref.obj,
    ref.ctype = 'CellType',
    marker.lst = marker.lst,
    method = 'SCINA',
    select.marker = 'Seurat',
    top.k = 30,
    strategy = 'marker-based'
) 
```

# Contributors
scAnnoX was developed by Xiaoqian Huang.