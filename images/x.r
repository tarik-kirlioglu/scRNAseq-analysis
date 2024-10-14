dn.markers <- c("S100A12", "S100A8", "S100A9", "PID1", "SASH1", "HLA-DMB", "CDKNL2",
                "FMNL2", "FCGR3A", "CLEC10A", "FCER1A", "ENHO", "CD3D", "CD3E", "CD3G",
                "CD8A", "CD8B", "TRGC2", "KLRF1", "CLIC3", "BNC2", "XCL1", "ZMAT4", 
                "NCAM1", "NELL2", "LEF1-ASI", "LEF1", "TDHZ2", "CCR7", "INPP4B", "FAAH2",
                "ICOS", "RTKN2", "IL2RA", "CTLA4", "GZMK", "IL7R", "KLRB1", "IGHD",
                "COL19A1", "TCL1A", "TNFRSF13B", "SOX5", "SSPN", "LILRA4", "SCT", "LRRC26",
                "TNFRSF17", "IGHA1", "IGLC2", "GP9", "CMTM5", "TMEM40")

DotPlot(object = seurat_obj, features = dn.markers) + RotatedAxis() + scale_color_gradient(low = "blue", high = "red")


markers <- c("PF4", "PPBP", #mega
             "CD79B", "CD79A", "MS4A1", "CD19", "CD72", "CD40", "CD22", "IGHD", "COL19A1", "TCL1A", "TNFRSF13B", "SOX5", "SSPN", #b cell
             "LDHB", "CD8B", "CD8A", "CD3D", "CD3E", "CD3G", "CD27", "IL7R", "CCR7", "NELL2", "LEF1", "ICOS", "RTKN2", "TSHZ2", "INP4B", "IL2RA", # t cell
             "LILRA4", "SCT", "LRRC26", "CLIC3", #pdc
             "CLEC10A", "FCER1A", "ENHO", #mdc
             "NKG7", "GZMB", "FGFBP2", "SPON2", "KLRC1", "FCGR3A", "GNLY", "KLRF1", "BNC2", "XCL1", "TRGC2", "ZMAT4", "NCAM11",  #nk cell
             "MS4A7", "FCN1", "S100A12", "S100A8", "CD14", "LYZ", "PID1", "FCGR3A", "FMNL2", #mono
             "GP9", "CMTM5", "TMEM40", #Platets
             "TNFRSF17", "TNFRSF13B", "IGHA1", "IGLC2", "SSPN") #plasmacells 