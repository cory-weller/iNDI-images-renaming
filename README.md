#!/usr/bin/env bash



```R
# after module load R/4.3

library(xlsx)
library(data.table)

map_filename <- '20240724_Renaming_image_map.xlsx'
if(!file.exists('formatted_map_60X.tsv')) {
    dat <- xlsx::read.xlsx(map_filename, sheetIndex=3, header=TRUE)
    dat.sheet1 <- xlsx::read.xlsx(map_filename, sheetIndex=1, header=TRUE)
    dat <- xlsx::read.xlsx(map_filename, sheetIndex=3, header=TRUE)
    setDT(dat)
    setnames(dat, c('Magnification','PlateID','WellID','Celltype','Batch','CellLine','Panel','Condition','Repetition','Markers','SiteIndex'))
    # 405, 488, 561, 640 nm

    dat[, paste0('well', 1:8) := tstrsplit(WellID, split=',')]
    dat[, paste0('Repetition', 1:8) := tstrsplit(Repetition, split=', ')]
    dat[, paste0('Marker_', c('405','488','561','640')) := tstrsplit(Markers, split=', ')]

    dat2 <- melt(dat, measure.vars=paste0('Repetition', 1:8), variable.name='Rep')
    dat2[, value := NULL]
    dat3 <- melt(dat2, measure.vars=paste0('well', 1:8), value.name='Well')
    dat3[, variable := NULL]
    dat3[, Repetition := NULL]
    dat3[, WellID := NULL]
    dat3[, Markers := NULL]
    dat3[, paste0('SiteIndex_', 1:9) := tstrsplit(SiteIndex, split=', ')]
    dat4 <- melt(dat3, measure.vars=paste0('SiteIndex_', 1:9), value.name='Site')
    dat4[, variable := NULL]
    dat4[,SiteIndex := NULL]
    setnames(dat4, 'Site','SiteIndex')

    # capitalize 'X' in magnification, and append 'air' to 20X
    dat4[Magnification=='20x', Magnification := '20Xair']
    dat4[Magnification=='60x', Magnification := '60X']

    # remove 'batch' from batch name, leaving only the integer
    dat4[, Batch := gsub('batch','',Batch)]

    # Capitalize 'P' in panel
    dat4[, Panel := gsub('panel','Panel',Panel)]

    # Capitalize 'D' at end of PlateID
    dat4[, PlateID := gsub('d$','D', PlateID)]
    dat4[, Column := gsub('[A-z]','',Well)]
    dat4[, Column := as.numeric(Column)]

    # Fix rows that should have batch 2 label
    dat4[CellLine %like% '2', Batch := 2]

    dat4 <- dat4[Magnification == '60X']
    dat4[, SiteIndex := NULL]
    dat4 <- unique(dat4)
    dat4[, idx := .I]
    dat4 <- merge(dat4, CJ('idx'=1:nrow(dat4), 'SiteIndex'=paste0('s',1:25)), by='idx')
    dat4[, idx := NULL]
    dat4[, Rep := NULL]
    dat4 <- unique(dat4)

    topdir <- '/data/INDIimage/Organellomics_Project/FUS-Pilot/FUS_D8_batch1-3/'
    dat4[, multichannel_fn := paste0(topdir, Magnification, '_OME-TIFFs/D8_batch', Batch, '/', 
                                    gsub('air','',Magnification), 'D8iNeuron_B', Batch, '_', Panel, '_', PlateID, '_', Well, '_', SiteIndex, '.ome.tif')]

    dat <- copy(dat4)
    rm(dat2)
    rm(dat3)
    rm(dat4)


    #fwrite(dat4, file='formatted_map_60X.tsv', quote=F, row.names=F, col.names=T, sep='\t')

# dat <- fread('formatted_map_60X.tsv')


# Cell types (e.g., U2OS cells / cortical neurons /microglia).
#     Batch: an experiment ID.
#         Cell lines (e.g., KOLF, TDP43, FUS/OPTN) 
#             Panel an panel ID (the 3-4 antibodies that were imaged together in the well)
#                 Conditions: (e.g., stress 1/ stress 2 / untreated )
#                     Reps (well in the plate). 
#                 Marker1
#                 Data (100 tiff images)
#                 Marker2
#                 Data (100 tiff images)
#                 Marker3
#                 Data (100 tiff images)
#                 DAPI
#                 Data (400 tiff images)
dat[, Celltype := NULL]
dat[, CELLTYPE := 'day8neuron']
dat[, BATCH := paste0('batch', Batch)]
setnames(dat, 'CellLine', 'CELLLINE')
dat[, CELLLINE := gsub('^FUS-H', 'FUSH', CELLLINE)]
dat[, CELLLINE := gsub('^FUS-R', 'FUSR', CELLLINE)]
dat[, CELLLINE := gsub('^FUS2-H', 'FUS2H', CELLLINE)]
dat[, CELLLINE := gsub('^FUS2-R', 'FUS2R', CELLLINE)]
setnames(dat, 'Panel', 'PANEL')
setnames(dat, 'Condition', 'CONDITION')

# The panel folders' names should start with a small letter, i.e. panelA, panelB,... instead of PanelA, PanelB
dat[, PANEL := gsub('Panel','panel',PANEL)]



# Capitalize
dat[, Marker_405 := gsub('Dapi','DAPI', Marker_405)]

dat <- melt(dat, measure.vars=c('Marker_405','Marker_488','Marker_561','Marker_640'), value.name='MARKER', variable.name='CHANNEL')
dat[CHANNEL=='Marker_405', CHANNEL := '0']
dat[CHANNEL=='Marker_488', CHANNEL := '1']
dat[CHANNEL=='Marker_561', CHANNEL := '2']
dat[CHANNEL=='Marker_640', CHANNEL := '3']
dat[, Rep := rleid(Well), by=list(CELLTYPE, BATCH, CELLLINE, PANEL, CONDITION, MARKER)]
dat[, Rep := as.character(Rep)]
# Rep names should be lower-case rep1-rep8
dat[, Rep := paste0('rep', Rep)]

dat[, singlechannel_fn := paste(CELLTYPE, BATCH, CELLLINE, PANEL, CONDITION, Rep, MARKER, sep='/')]
dat[, singlechannel_fn := paste0(singlechannel_fn, '/', SiteIndex, '.tif')]

fwrite(dat, file='formatted_60X_map.csv', quote=F, row.names=F, col.names=T, sep=',')



```
```bash
cp /data/INDIimage/Organellomics_Project/FUS-Pilot/FUS_D8_batch1-3/20Xair_OME-TIFFs/D8_batch1/20XD8iNeuron_B1_PanelA_indi000000102D_C3_s1.ome.tif .

wget https://downloads.openmicroscopy.org/bio-formats/6.0.1/artifacts/bftools.zip
unzip bftools.zip
alias bfconvert='/vf/users/CARD/users/wellerca/iNDI-images-renaming/bftools/bfconvert'

bfconvert -channel 0 image.tif image.0.tif
bfconvert -channel 1 image.tif image.1.tif
bfconvert -channel 2 image.tif image.2.tif
bfconvert -channel 3 image.tif image.3.tif

D28nuron_batch2-p15_FUSRev_PanelA_untreated_C10Pt-1_Dapi(405),FMRP(488),SQSTM1(561),TUJ1(640)


`
- in `/data/INDIimage/`, rename `Organellomics Project`  -> `Organellomics_Project`

# find . -name '*.tif' -exec rename '_B1-2_' '_B2_' {} \;