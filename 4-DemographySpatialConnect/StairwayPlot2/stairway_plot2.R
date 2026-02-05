library(dartR.popgen)

setwd('~/Documents/guzmanrs/Paxil_proj/stairway_plot/')

paxil_gl <- gl.read.vcf('../../1-VariantCallFilt/07_freebayes/Paxil_M095_noLD.recode.vcf.gz', mode='genotype')

paxil_info <- read.csv('../../Paxil_PopGroupMap.csv')

row.names(paxil_info) <- paxil_info$Indv
paxil_info <- paxil_info[paxil_gl@ind.names,]

pop(paxil_gl) <- paxil_info$Pop

pops <- unique(paxil_info$Pop)

paxil_singlepops_gls <- list()
for (i in 1:length(pops)){
  paxil_singlepops_gls[[i]] <- gl.drop.pop(paxil_gl, pop.list = pops[-c(i)])
}
names(paxil_singlepops_gls) <- pops

sw_grps <- list()
st.path <- 'stairway_plot_v2.1.1/'
plotdir <- './'
for (i in 1:length(pops)){
  sw_grps[[i]] <- gl.run.stairway2(
    paxil_singlepops_gls[[i]], 
    L=1218700, 
    mu=1e-8, 
    nreps=1000, 
    parallel=14, 
    stairway2.path = st.path, 
    stairway_plot_dir = paste0(st.path,"/stairway_plot_es"),
    plot.dir = plotdir, 
    filename = ''
  )
}

rds_files <- list.files(path = '.', pattern = '.RDS')
svg_files <- list.files(path = '.', pattern = '.svg')
rds_cont <- list()

svg('paxil_stairwayplot_051224_Axil03.svg', width = 8, height = 8)
readRDS('paxil_nomiss_051224_Axil03.RDS')
dev.off()

svg('paxil_stairwayplot_051224_Axil04.svg', width = 8, height = 8)
readRDS('paxil_nomiss_051224_Axil04.RDS')
dev.off()

svg('paxil_stairwayplot_051224_Axil05.svg', width = 8, height = 8)
readRDS('paxil_nomiss_051224_Axil05.RDS')
dev.off()

svg('paxil_stairwayplot_051224_Axil06.svg', width = 8, height = 8)
readRDS('paxil_nomiss_051224_Axil06.RDS')
dev.off()

svg('paxil_stairwayplot_051224_Axil07.svg', width = 8, height = 8)
readRDS('paxil_nomiss_051224_Axil07.RDS')
dev.off()

svg('paxil_stairwayplot_051224_Axil08.svg', width = 8, height = 8)
readRDS('paxil_nomiss_051224_Axil08.RDS')
dev.off()

svg('paxil_stairwayplot_051224_Axil09.svg', width = 8, height = 8)
readRDS('paxil_nomiss_051224_Axil09.RDS')
dev.off()

svg('paxil_stairwayplot_051224_Axil10.svg', width = 8, height = 8)
readRDS('paxil_nomiss_051224_Axil10.RDS')
dev.off()

svg('paxil_stairwayplot_051224_Axil11.svg', width = 8, height = 8)
readRDS('paxil_nomiss_051224_Axil11.RDS')
dev.off()

svg('paxil_stairwayplot_051224_Axil43.svg', width = 8, height = 8)
readRDS('paxil_nomiss_051224_Axil43.RDS')
dev.off()

svg('paxil_stairwayplot_051224_Axil63.svg', width = 8, height = 8)
readRDS('paxil_nomiss_051224_Axil63.RDS')
dev.off()

svg('paxil_stairwayplot_051224_AxilIS.svg', width = 8, height = 8)
readRDS('paxil_nomiss_051224_AxilIS.RDS')
dev.off()

svg('paxil_stairwayplot_051224_allPops.svg', width = 8, height = 8)
readRDS('paxil_nomiss_stairway_051224.RDS')
dev.off()

###### Groups
svg('paxil_stairwayplot_061224_Axil02.svg', width = 8, height = 8)
readRDS('paxil_nomiss_061224_Axil02.RDS')
dev.off()

svg('paxil_stairwayplot_061224_Center.svg', width = 8, height = 8)
readRDS('paxil_nomiss_061224_Center.RDS')
dev.off()

svg('paxil_stairwayplot_061224_CenterNorth.svg', width = 8, height = 8)
readRDS('paxil_nomiss_061224_CenterNorth.RDS')
dev.off()

svg('paxil_stairwayplot_061224_South.svg', width = 8, height = 8)
readRDS('paxil_nomiss_061224_South.RDS')
dev.off()

svg('paxil_stairwayplot_061224_North.svg', width = 8, height = 8)
readRDS('paxil_nomiss_061224_North.RDS')
dev.off()


