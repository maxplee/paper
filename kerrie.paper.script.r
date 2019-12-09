# R scripts for Kerrie L. Marie's paper "Melanoblast transcriptome analysis 
# reveals novel pathways promoting melanoma metastasis"
# Nature communication

Dir = '/data/LeeGrp/Yang/collab/glenn/Kerrie/share/'
f1 = 'E15.5E17.5P1P7_gns467_z'; cex1 = 0.2
E.P.M3M4.heatmaps(Dir, f1, cex1, use.log2=F, delta=0.6)
# E15.5E17.5P1P7_gns467_z_heatmap.pdf  for Fig.1c 
f1 = 'E15.5E17.5P1P7_gns46_z'; cex1 = 0.6
E.P.M3M4.heatmaps(Dir, f1, cex1, use.log2=F, delta=0.6)
# E15.5E17.5P1P7_gns46_z_heatmap.pdf for Fig.1d

Dir = '/data/LeeGrp/Yang/collab/glenn/Kerrie/share/GSE8401/'
f1 = 'GSE8401_log2RMA_gene'
f.w = 'E.vs.P.padj0.1GSE19234FC1.5P0.1HR.gt1.wt1' # 43 genes, p=3.5e-05
f.clinic = 'GSE8401surv_OS_month'
testGSE8401.wt.score.KM(Dir, f1, f.w, f.clinic, wt.name=f.w, cex=1)
# GSE8401_log2RMA_gene_E.vs.P.padj0.1GSE19234FC1.5P0.1HR.gt1.wt1_score_surv_Late.event.os.score.bin.KM.pdf  for Fig.1e 
# GSE8401_log2RMA_gene_E.vs.P.padj0.1GSE19234FC1.5P0.1HR.gt1.wt1_score_surv_Early.event.os.score.bin.KM.pdf for Fig.1f

Dir = '/data/LeeGrp/Yang/collab/glenn/Kerrie/share/GSE8401/'
f1 = 'E_P_GSE8401_log2RMA_gene79gns'; cex0=0.4
f2 = 'infSamples'
GSE8401.hm(Dir, f1, f2, cex0=cex0, mode='z.stats')
# E_P_GSE8401_log2RMA_gene79gns.z.stats.hm.pdf for Supp Figure 3a 

Dir = '/data/LeeGrp/Yang/collab/glenn/Kerrie/share/GSE98394/'

f1 = 'GSE98394_expression66out79_z.row2groups_wt_score_stage'
inName = paste(Dir, f1, sep='')
X = my.read.tb(inName)
z = split(X[,'score'], as.vector(X[,'class']))
tmp = kruskal.test(score ~ class, data=X)
pv = signif(tmp$p.value, digit=2)

grp.names = names(z); grp.names[1] = 'Nevus'
ylim = range(z); ylim[2] = 80
y.pv = c(45, 65, 75)
outName = paste(Dir, f1, '.boxplot.pdf', sep='')
pdf(file=outName)
run.beeswarm.boxplot(z, grp.names, ylim, y.pv, xlab='Stage', ylab='Score',
pv=pv, pv.desc='Kruskal-Wallis ')
dev.off()
# out GSE98394_expression66out79_z.row2groups_wt_score_stage.boxplot.pdf for
# Supp Figure 3b


E.P.M3M4.heatmaps =
function(Dir, f1, cex1, use.log2=T, upper=NULL, delta=1)
{

f2 = 'E15.5E17.5P1P7_gns46'
gn.vec = c('Kdelr3', 'Dab2', 'P4ha2', 'Gulp1')
inName = paste(Dir, f1, sep='')
X = my.read.tb(inName)
inName = paste(Dir, f2, sep='')
Y = my.read.tb(inName)
x = as.vector(X[,1])
map.ylab.col = rep('black', length(x))
s1 = is.element(x, as.vector(Y[,1]))
map.ylab.col[s1] = 'green'
s1 = is.element(x,  gn.vec)
map.ylab.col[s1] = 'magenta'
names(map.ylab.col) = x

inName = paste(Dir, f1, sep='')
X = my.read.tb(inName)
gns = as.vector(X[,1])
X = data.matrix(X[,-1])
rownames(X) = gns



if (use.log2)
Y = log2(1+X)
else
Y = X

if (!is.null(upper))
{
Y = my.saturation.matrix(Y, lower=0, upper)
outName = paste(Dir, f1,'_upper', upper,  '_heatmap.pdf', sep='')
}
else
outName = paste(Dir, f1, '_heatmap.pdf', sep='')

pdf(file=outName)
matrix.image(Y, ylabs=gns, cex0=0.9, row.hclust=T, put.tick=F,
colors='bwr', delta=delta,ylab.cex=cex1, line1=0.5, ylab2col=map.ylab.col)
dev.off()

cat('write', outName, '\n') 

}


GSE8401.hm =
function(Dir, f1, f2, cex0=0.4, mode='scale', col.hclust=F, at.vec=NULL)
{
# mode='scale' or mode='z.stats'

inName = paste(Dir, f1, sep='')
X = my.read.tb(inName)
gns = as.vector(X[,1])
X = data.matrix(X[,-1])
rownames(X) = gns
if (mode=='z.stats')
{
X[,1:4] = t(apply(X[,1:4], 1, z.stats))
X[,-(1:4)] = t(apply(X[,-(1:4)], 1, z.stats))
X = my.saturation.matrix(X, lower=-2, upper=2)
}
else
{
m1 = mean(as.vector(X[,1:4]))
v1 = var(as.vector(X[,1:4]))
X[,-(1:4)] = scale.data.matrix(X[,-(1:4)], m1, v1)
}

n = dim(X)[2]
f2 = 'infSamples'
y = GSE8401PT.MT.ids(Dir, f2)
i1 = y$PT
i2 = y$MT
cols = c(colnames(X[,1:4]), i1, i2)
X = X[,cols]
col.grp.num = c(rep(1,2), rep(2,2), rep(3, length(i1)), rep(4, length(i2)))
names(col.grp.num) = colnames(X)

col.vec = c('blue', 'green', 'yellow', 'red')
grp.names = c('E', 'P', 'Primary', 'Metastatic')
if (is.null(at.vec))
at.vec = c(0.5, 3, 0.2*n, 0.7*n)/n

if (col.hclust)
px1 = 'col.hclust.'
else
px1 =NULL

if (mode=='z.stats')
outName = paste(Dir, px1,  f1, '.', mode, '.hm.pdf', sep='')
else
outName = paste(Dir, px1, f1, '.hm.pdf', sep='')

pdf(file=outName)
matrix.image.col.grp.bar(X, col.grp.num, col.vec=col.vec,  cex0=cex0,
heights = c(9, 1), grp.names=grp.names, at.vec=at.vec, col.hclust=col.hclust)
dev.off()

}


testGSE8401.wt.score.KM =
function(Dir, f1, f.w, f.clinic, wt.name, sx1 = '_score', cex=1,
col.vec=c('blue', 'red'), ex4PTs=T)
{

z = run.weighted.index(Dir, f1, f.w)
X = data.frame(sample.id=names(z), score=z)
f2 = paste(f1, '_',  f.w, sx1, sep='')
outName = paste(Dir, f2, sep='') 
write.table(X, file=outName, sep='\t', quote = F, row.names=F)
cat('write', outName, '\n') 

inName = paste(Dir, f.clinic, sep='')
X = my.read.tb(inName)
 table(X[,'SampleType'], X[,'Stage2'])
# to exclude 4 primary samples in Late stage.
s1 = (X[,'SampleType']=='Primary Site')&(X[,'Stage2']=='Late')
if (ex4PTs)
X = X[!s1,]

inName = paste(Dir, f2, sep='')
Y = my.read.tb(inName)
X = merge(X, Y, by='sample.id')
f3 = paste(f2, '_surv', sep='')
outName = paste(Dir, f3, sep='')
write.table(X, file=outName, sep='\t', quote = F, row.names=F)
cat('write', outName, '\n') 

stages = c('Early', 'Late')
for (v1 in stages)
{
Y = subTable(X, 'Stage2', v1)
outName = paste(Dir, f3, '_', v1, sep='')
write.table(Y, file=outName, sep='\t', quote = F, row.names=F)
cat('write', outName, '\n') 
}

gs = paste(f3, '_', stages, sep='')
gs = c(gs, f3)
for (g1 in gs)
{
inName = paste(Dir, g1, sep='')
Y = my.read.tb(inName)
y = cut.bin.vrb(Y[,'score'], prob=0.5,  mode='0:1')
Y = cbind(Y, score.bin=y)
event.time = 'time.os'; event = 'event.os'; ylab='Overall Survival'
vrb1 = 'score.bin'
cols = c(event.time, event, vrb1)
grp.names = c('Low', 'High')
names(grp.names) = c('Group0', 'Group1')
if (regexpr('Early', g1)>0)
x1 = 50
else
x1 = 0

outName = paste(Dir, g1, '.', event,'.', vrb1, '.KM.pdf',  sep='')
pdf(file=outName)
my.survfit(Y, cols, col.vec=col.vec, xy=c(x1+30,0.9), p.xy=c(x1+38, 0.7),
xlab='Month', ylab=ylab, grp.names=grp.names, cex=cex)
dev.off()

}

}


run.beeswarm.boxplot =
function(y, grp.names, ylim, y.pv=NULL, ylab=NULL, xlab=NULL, pv=NULL,
pv.desc=NULL, col.vec=c('black', 'green', 'blue', 'magenta'), put.pv=T,
method='swarm', cex.pt=0.8)
{
# y.pv, y coordinates to put the p-values.
# Note: Before running this fun, do boxplot(y) to determine y coordinates
# to put the p-values and the grp.names for the boxplots. 
# The Wilcoxon test was used to compare paired groups to compute p-values.

beeswarm.boxplot(y, pv.desc=pv.desc, pv=pv, xlab=xlab, ylab=ylab,
ylim=ylim, names=grp.names, col.vec=col.vec, cex.axis=1.3, cex.lab=1.3, 
cex.pt=cex.pt, put.pv=put.pv, method=method)

pv.vec = signif(wilcox.test2levels(y), digit=2)

n = length(grp.names)
a = cbind(1:(n-1), 2:n)

if (!is.null(y.pv))
{
y = y.pv
for (i in 1:(n-1))
{
my.bridge(a[i,], y[i])
text(a[i,1]+0.5, y[i]+4, pv.vec[i], cex=1.2)
}
}

}


my.saturation.matrix =
# input is a data matrix without any missing values
function(X, lower, upper)
{
S1 = X<lower
X[S1] = lower
S1 = X>upper
X[S1] = upper

return(X)
}


matrix.image =
function(X, cex0=0.8, x.las=2, ylabs=NULL, left=4, put.tick=T, ncolor=100,
row.hclust=F, colors='topo', zlim=NULL, put.xlab=T, col.hclust=F, delta=2,
xlab.cex=NULL, ylab.cex=NULL,  line1=0.5, bottom=5, top=4, right=5, xlab2col=NULL,
add.grid=F, ylab2col=NULL, bias=1)
{
# set zlim so that white corresponding to 0: zlim = c(-v1, v1)
# or set zlim for multiple heatmaps with same color bar legends.
# Use delta to position the image bar
# use colors='bwr.bias' and bias=0.5 or 2 for the data with an asymetric 
# distribution.


if (is.null(xlab.cex))
xlab.cex = cex0

if (is.null(ylab.cex))
ylab.cex = cex0


if (row.hclust)
{
hc2 = hclust(dist(X))
rowOrder = hc2$order
X = X[rowOrder,]
if (!is.null(ylabs))
ylabs = ylabs[rowOrder]
}

if (col.hclust)
{
hc = hclust(dist(t(X)))
colOrder = hc$order
X = X[,colOrder]
}

if (colors=='topo')
col.vec = topo.colors(ncolor)

if (colors=='bwr')	# blue-white-red colors
col.vec = colorRampPalette(c("blue", "white", "red"), space = "rgb")(101)

if (colors=='bwr.bias')  # blue-white-red colors
col.vec = colorRampPalette(c("blue", "white", "red"), space = "rgb", bias=bias)(101)


if (colors=='wr')  # white-red colors
col.vec = colorRampPalette(c("white", "red"), space = "rgb")(101)

if (colors=='wb')  # white-blue colors
col.vec = colorRampPalette(c("white", "blue"), space = "rgb")(101)

if (colors=='gyr')  	# green-yellow-red colors
col.vec = colorRampPalette(c("green", "yellow", "red"), space = "rgb")(101)
if (colors=='bgr')	# blue-gray-red colors
col.vec = colorRampPalette(c("blue", "gray", "red"), space = "rgb")(101)
if (colors=='gkr')		# green-black-red colors
col.vec = colorRampPalette(c("green", "black", "red"), space = "rgb")(101)

n = ncol(X)
x = 1:n
m = nrow(X)
y = 1:m
if (is.null(zlim))
zlim = range(X, na.rm=T)	 #zlim = signif(range(X, na.rm=T), digit=2)
par(mar=c(bottom, left,top, right), xpd=NA)
image(x,y,t(X), col=col.vec, axes = FALSE, ylab='', xlab='')
if (add.grid)
my.add.grid(x,y)

image.legend(n+delta, m, zlim, col=col.vec, cex=0.5)
xlabs = colnames(X)
if (put.xlab)
mtext(xlabs, side=1,at=x,line=0.5,cex=xlab.cex, las=x.las, col=xlab2col[xlabs])
#axis(1, labels=xlabs, at=x, line=line1, cex.axis=xlab.cex, las=x.las, tick=put.tick, col=xlab2col[xlabs] )
# change line1 to change the label position


if (!is.null(ylabs))
#axis(2, labels=ylabs, at=y, line=line1, cex.axis=ylab.cex, las=1, tick=put.tick, col.axis = col.ylab)
mtext(ylabs, side=2,at=y,line=line1, cex=ylab.cex, las=1, col=ylab2col[ylabs])

return(X)
}


my.read.tb =
function(inName, nrows=-1, quote = "\"", comment.char = "#", header=T, 
check.names=F, na.strings="NA", skip=0, sep='\t')
{
# check.names is F and sep is '\t'
X = read.table(inName, row.names=NULL, header=header, sep=sep,
check.names=check.names, nrows=nrows, quote=quote, comment.char=comment.char,
na.strings=na.strings, skip=skip)
}


GSE8401PT.MT.ids =
function(Dir, f1)
{

inName = paste(Dir, f1, sep='')
Y = my.read.tb(inName)
s1 = !is.na(Y[,'Stage2'])
Y = Y[s1,]
s1 = (Y[,'SampleType']=='Metastasis')&(Y[,'Stage2']=='Late')
i2 = as.vector(Y[s1,1])
s1 = (Y[,'SampleType']=='Primary Site')&(Y[,'Stage2']=='Early')
i1 = as.vector(Y[s1,1])

y = list(PT=i1, MT=i2)

return(y)

}


matrix.image.col.grp.bar =
function(X, col.grp.num, col.vec = c('green', 'magenta'), cex0=0.6,
heights = c(8, 2), grp.names=c('Nevus', 'primary melanoma'),
at.vec=c(13, 52)/78, ylab2col=NULL, col.hclust=F, delta=2, ylab.cex=0.5)
{
# This fun add the bars under the heatmap to indicate the groups.
# col.grp.num is a named vector to map column names in X to a vector of
# integers representing groups.
# delta is a parameter to change the position of the image bar.

def.par <- par(no.readonly = TRUE) # save default, for resetting...

layout(mat=matrix(1:2,2,1,byrow=TRUE),  heights=heights, TRUE)

gns = rownames(X)
X = matrix.image(X, ylabs=gns, cex0=cex0, row.hclust=T,
put.tick=F, colors='bwr', delta=delta, ylab.cex=ylab.cex,
put.xlab=F, col.hclust=col.hclust, bottom=0, top=3, right=5, ylab2col=ylab2col)

y = col.grp.num[colnames(X)]

par(mar=c(2,4,0.5,5))
Y = data.matrix(y)
image(Y, col=col.vec, axes = FALSE)

if (col.hclust)
{
mtext(grp.names, side = 1, col=col.vec, at =at.vec)
}
else
mtext(grp.names, side = 1, at = at.vec)

par(def.par)  #- reset to default

}


run.weighted.index =
function(Dir, f1, f.w, ex.cols=1)
{
inName = paste(Dir, f1, sep='')
X = my.read.tb(inName)

if (is.character(f.w))
{
inName = paste(Dir, f.w, sep='')
Y = my.read.tb(inName)
}
else
Y = f.w

s1 = is.element(as.vector(X[,1]), as.vector(Y[,1]))
X = X[s1,]
gns = as.vector(X[,1])
X = data.matrix(X[,-ex.cols])
rownames(X) = gns
z = weighted.index(X, Y)

return(z)

}

my.survfit =
function(X, cols, col.vec=c('red',  'blue', 'magenta', 'cyan'),
title1=NULL, xy=c(2,0.25), p.xy=c(5, 0.4), ylab='survival probability',
xlab='year', cex=0.8, yscale=1, KM.plot=T, grp.lab='Group',
my.mar=c(5, 4.5, 4, 2), add.HR=F, grp.names=NULL, lwd=1, solid.line=F)
{

# 4 colors in col.vec are good for 4 or less number of values in X[,vrb].
# grp.lab='Q' for groups Q1, Q2, Q3 and Q4.
# grp.names = c('Low index', 'High index')
# names(grp.names) = c('Group0', 'Group1')

library(survival)

names(cols) = c('time', 'cens', 'variable')
s1 = !is.na(X[,cols['time']])
X = X[s1,]

vrb = cols['variable']
x = X[,vrb]
s1 = (x=='')|(x=='  ')|(x=='Unclassified')
x[s1] = NA
nv = lenu(excludeNA(x))		# number of non-NA values
col.vec = col.vec[1:nv]

y = Surv(X[,cols['time']], X[,cols['cens']])
surv1 = survfit(y ~ x)
diff1 = survdiff(y ~ x)
chisq1 = diff1$chisq
df1 = length(diff1$n)-1
pv = 1 - pchisq(chisq1, df=df1)
# pv = min(1, 2*pv)		# to be consistent with the log-rank pv in my.logrank.test().
pv = signif(pv, digit=4)
if (add.HR)
{
coxph1 = coxph(y~x)
tmp = coxph.stats(coxph1, num.signif=4)
w = hr.CI.coxph1variable(tmp, cols[3], digit=2)
print(w)
}

if (solid.line)
lty.vec = rep(1, nv)
else
lty.vec = 1:nv

if (KM.plot)
{
par(mar=my.mar)
plot(surv1, lty=lty.vec,  col=col.vec, xlab=xlab, ylab=ylab, yscale=yscale,
cex=cex,cex.axis=cex, cex.lab=cex, mark.time=T, lwd=lwd)
title(title1)
if (add.HR)
str2 = paste('p-value=', signif(pv, digit=4),  ' HR=', w['HR'], sep='') 
else
str2 = paste('p-value=', pv, sep='') 

sizes = table(x)
x1 = names(surv1$strata)
x1 = sub('x=', '', x1)
x2 = sizes[x1]
if (grp.lab=='Q')
{
map2grp = c('Q1', 'Q2', 'Q3', 'Q4')
names(map2grp) = as.character(c(-2, -1, 1, 2))
x1 = map2grp[x1]
}
else
x1 = paste(grp.lab, x1, sep='') 

if (!is.null(grp.names))
x1 = grp.names[x1]

legend.lab = paste(x1, x2, sep=' n=')
legend(xy[1], xy[2], legend=legend.lab, lty=lty.vec, col=col.vec,
text.col=col.vec,cex=cex)
text(p.xy[1], p.xy[2], str2, cex=cex)
}

if (add.HR)
{
z = c(w, pv)
names(z) = c('gene', 'HR', 'CI95percent', 'Coxph.pv', 'KM.pv')
return(z)
}
else
return(pv)

}

beeswarm.boxplot =
function(y, ylab=NULL, xlab=NULL, extend1=NULL, lwd=1.5, pv=NULL, pv.desc=NULL,
cex.axis=1, cex.lab=1,  col.vec=c('skyblue', 'purple'), cex.pt=0.8, put.pv=T,
use.beeswarm = T, put.size=T, ylim=NULL, names=NULL, method='swarm')
{

library(beeswarm, lib.loc='/home/yanghow/R/x86_64-unknown-linux-gnu-library/3.0/')

if (is.null(ylim))
ylim = range(y)

if (is.null(names))
names = names(y)

if (use.beeswarm)
{
boxplot(y, ylab=ylab, xlab=xlab, lwd=lwd, cex.axis=cex.axis, cex.lab=cex.lab,
ylim=ylim, names=names)
beeswarm(y, add=T, col=col.vec, cex=cex.pt, method=method, corral='gutter')
}
else
boxplot(y, ylab=ylab, xlab=xlab, lwd=lwd, cex.axis=cex.axis, cex.lab=cex.lab, 
col=col.vec)


if (is.null(pv))
{
tmp = t.test(y[[1]], y[[2]])
pv = tmp$p.value
}

if (put.pv)
{
m1 = paste(pv.desc, ' p=', signif(pv, digit=4), ' ', extend1, sep='')
mtext(m1, cex=cex.lab)
}

if (put.size)
{
m = length(y)
n.vec = c()
for (i in 1:m)
n.vec = c(n.vec, length.ex.na(y[[i]]))
str.vec = paste('n=', n.vec, sep='')
mtext(str.vec, side=1, line=2, at=1:m, cex=cex.lab)

}

return(pv)

}
