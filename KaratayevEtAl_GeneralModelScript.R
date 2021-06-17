require(mgcv); require(fields);
ddall=read.csv("DreissenaGreatLakes1990_2018_Pub.csv")
ddall$Lake=as.character(ddall$Lake)
ddall$Lake[ddall$Lake=="Erie" & ddall$Basin=="Western"]="ErieWest";
ddall$Lake=as.factor(ddall$Lake)
colnames(ddall)[10:15]=c("ZMN","QMN","DN","ZMB","QMB","DB")



#Model comparison function
ModelComp=function(Lk,type=2,MaxKs=c(16,25)){
	Dat=na.omit(ddall[ddall$Lake%in%Lk,c(1,2,5,15)])
	FAM=tw(); if(type==1){ FAM=gaussian(); Dat$DB=log(1+Dat$DB); };
	ALL=gam(DB~s(Year,Depth.m,k=MaxKs[2])+as.factor(Lake),data=Dat,family=FAM)
	LYD=gam(DB~s(Year,k=MaxKs[1])+s(Depth.m,k=MaxKs[2])+as.factor(Lake),data=Dat,family=FAM)	
	YD=gam(DB~s(Year,k=MaxKs[1])+s(Depth.m,k=MaxKs[2]),data=Dat,family=FAM)
	YDI=gam(DB~s(Year,Depth.m,k=MaxKs[2]),data=Dat,family=FAM)
	LY=gam(DB~s(Year)+as.factor(Lake),data=Dat,family=FAM)
	LD=gam(DB~s(Depth.m,k=MaxKs[2])+as.factor(Lake),data=Dat,family=FAM)
	D=gam(DB~s(Depth.m,k=MaxKs[2]),data=Dat,family=FAM)
	Y=gam(DB~s(Year,k=MaxKs[1]),data=Dat,family=FAM)
	L=gam(DB~as.factor(Lake),data=Dat,family=FAM)
	rsq=function(x) cor(log(1+Dat$DB),log(1+x$fitted.values))^2
	ResSum=round(cbind(round(AIC(ALL,LYD,LY,LD,YD,YDI,D,Y,L)),"R-squared"=c(rsq(ALL),rsq(LYD),rsq(LY),rsq(LD),rsq(YD),rsq(YDI),rsq(D),rsq(Y),rsq(L))),2)
	Final=ResSum[order(ResSum[,2]),c(2,1,3)]; return(Final);
}


#Plotting functions:
zCatPlot=function() abline(h=c(30,50,90),col="white",lwd=3,lty=3)
DomPlot=function(Lk){
	dtmp=ddall[ddall$Lake%in%Lk,]; asn=function(x) as.numeric(as.character(x)); naz=function(x){ x[is.na(x)]=0; return(x); }; 
	QD=as.matrix(aggregate(cbind(0,asn(dtmp$QMB),asn(dtmp$DB),asn(dtmp$ZMB)),by=list(dtmp$Year),FUN=sum,na.rm=TRUE)); QD[,2]=smooth(naz(QD[,3]/QD[,4]));
	plot(QD[,1:2],ylim=c(0,1),xlim=c(QD[1,1]+0.5,2018),type="l",main=paste0(Lk,collapse=" ")); 
	YrSq=c(QD[,1],rev(QD[,1])); polygon(YrSq,c(0*QD[,1],rev(QD[,2])),border=0,col="brown2"); polygon(YrSq,c(1+0*QD[,1],rev(QD[,2])),border=0,col="grey40");
}
clearR=10; Clear=function() pmin((ddall$DB*1e3*clearR*24)/(ddall$Depth.m*1e6),1)

#Smoothing surface function. k sets surface smoothness. MGCV defaults to k=25. 
#We use smaller k (see below), especially in shallow lakes with fewer data points.
smoo2d2=function (dat, k = 15, mlk=FALSE, cats = 1:3, n.grid = 30, nest = NULL, 
    rngMaxes = NULL, plotgive = FALSE, nonzeroCenter = TRUE, 
    sampPt = list(pch = 16, cex = 0.25, col = 1, jitter = 0)) 
{
    n.grid=c(head(n.grid,1),tail(n.grid,1))
    k=c(head(k,1),tail(k,2)[1],tail(k,1))
	datDF = data.frame(dat)
    datnm = names(datDF)
    sampts = function() points(jitter(dat[, cats[-1]], sampPt[[4]]), 
        pch = sampPt[[1]], cex = sampPt[[2]], col = sampPt[[3]])
	Form=paste(datnm[1], "~", "s(",datnm[cats[2]],",",datnm[cats[3]],",k=",k[3],")")
	if(mlk) Form=paste(Form,"+s(",datnm[cats[2]],",k=",k[1],")",  "+s(",datnm[cats[3]],",k=",k[2],")")
	if(length(cats)==4){
		PRM=as.factor(dat[,cats[4]]); parMns=as.matrix(aggregate(dat[,1]~PRM,FUN=function(x) c(mean(x,na.rm=TRUE),length(x)/nrow(dat)))); parMns=apply(parMns,2,function(x) as.numeric(as.character(x)));
		datDF[,1]=datDF[,1]-parMns[match(datDF[,4],parMns[,1]),2]; Gmean=sum(parMns[,2]*parMns[,3]);
	}	
	xgam=gam(as.formula(Form),data = datDF)

 if (plotgive == 1) {
        vis.gam(xgam, plot.type = "contour", color = "terrain", 
            contour.col = "grey40", lwd = 1, labcex = 0.1)
        sampts()
    }
    r1 = range(xgam$var.summary[[1]])
    r2 = range(xgam$var.summary[[2]])
    if (!is.null(rngMaxes)) {
        r1[2] = rngMaxes[1]
        r2[2] = rngMaxes[2]
    }
    newd = data.frame(matrix(0, prod(n.grid), 0))
    newd[[1]] = rep(seq(r1[1], r1[2], len = n.grid[1]), n.grid[2])
    newd[[2]] = rep(seq(r2[1], r2[2], len = n.grid[2]), rep(n.grid[1],n.grid[2]))
	names(newd) = names(xgam$var.summary);
    z = t(matrix(predict.gam(xgam, newdata = newd, type = "link"), n.grid[1], n.grid[2]))
	if(length(cats)>3) z=z+Gmean
    if (nonzeroCenter) 
        z = pmax(z, 0) * mean(dat[, 1])/mean(pmax(z, 0))
    rown = round(unique(newd[[2]]), 2)
    coln = round(unique(newd[[1]]), 2)
    if (!is.null(nest)) {
        if (nest == cats[2]) 
            z = z %*% diag(coln/colMeans(z))
        else z = t(t(z) %*% diag(rown/rowMeans(z)))
    }
    colnames(z) = coln
    rownames(z) = rown
    if (plotgive > 1) {
        SpatiotempPlot(pmax(z, 0), Xax = coln, Yax = rown, palette = plotgive - 1)
        sampts()
    }
    return(z)
}
#Function to plot heatmaps
SpatiotempPlot=function(spacetime, Xax = 1:dim(spacetime)[2], Yax = 1:dim(spacetime)[1], 
    XaxN = "X", YaxN = "Y", figtitle = "Title", Zlim = c(0, max(spacetime, 
        na.rm = TRUE)), cont = NULL, cexAx = 1, contPlot = spacetime, 
    cexCont = 1.5 * cexAx, Conts = NULL, contSpan = 1, palette = 1) 
{
    require(fields)
    spacetime[is.na(spacetime)] = Zlim[1] - 1
    COL = rev(rainbow(1000, start = 0, end = 0.7))
    if (palette > 1) {
        require(pals)
        COL = list(parula(1000), head(tail(parula(1000), -50, 
            -50)))[[palette - 1]]
    }
    if (length(Zlim) == 1) {
        Zlim = quantile(spacetime, c(Zlim, 1 - Zlim), na.rm = T)
        Rng = range(spacetime)
        if (Zlim[1] < Rng[1]) 
            Zlim[1] = Rng[1]
        if (Zlim[2] > Rng[2]) 
            Zlim[2] = Rng[2]
    }
    spacetime[which(is.na(spacetime), arr.ind = TRUE)] = max(Zlim) + 
        1
    image.plot(x = Xax, y = Yax, z = t(spacetime), zlim = Zlim, 
        xlab = XaxN, ylab = YaxN, cex.axis = cexAx, cex.lab = cexAx, 
        legend.cex = cexAx, main = figtitle, col = COL)
    box()
    if (!is.null(cont)) {
        if (abs(log(max(contPlot, na.rm = T), 10)) > 2) 
            options(scipen = -10)
        if (contSpan > 1) {
            smoo1 = t(apply(contPlot, 1, function(x) supsmu(1:ncol(contPlot), 
                x, span = contSpan/ncol(contPlot))$y))
            smoo2 = t(apply(smoo1, 1, function(x) supsmu(1:ncol(smoo1), 
                x, span = contSpan/ncol(smoo1))$y))
            contPlot = smoo2
        }
        if (is.null(Conts)) 
            contour(x = Xax, y = Yax, z = t(contPlot), add = T, 
                col = cont, lwd = 1.5 * cexCont, labcex = cexCont)
        if (!is.null(Conts)) 
            contour(x = Xax, y = Yax, z = t(contPlot), levels = Conts, 
                add = T, col = cont, lwd = 1.5 * cexCont, labcex = cexCont)
        if (abs(log(max(contPlot, na.rm = T), 10)) > 2) 
            options(scipen = 0)
    }
}


#Run analyses:
#Deep lakes
Lk=c("Huron","Ontario","Michigan"); 
#model comparison table
ModelComp(Lk,MaxKs=c(years=10,depths=10),type=2);
#Mussel biomass over depths and time
x=smoo2d2(k=25,nonzeroCenter=FALSE,rngMaxes=c(2018,200),n.grid=40,cats=1:3,plotgive=2,sampPt=list(pch=1,cex=0.5,col="darkgrey",jitter=0),dat=na.omit(cbind(log(1+ddall$DB),ddall$Year,ddall$Depth.m,ddall$Lake)[ddall$Lake%in%Lk,])); zCatPlot();
#Species composition
DomPlot(Lk);
#Mussel filtering intensity over depths and time
x=smoo2d2(k=25,nonzeroCenter=FALSE,rngMaxes=c(2018,200),n.grid=40,cats=1:3,plotgive=2,sampPt=c(pch=1,cex=0.5,col=1,jitter=0),dat=na.omit(cbind(Clear(),ddall$Year,ddall$Depth.m,ddall$Lake)[ddall$Lake%in%Lk,])); zCatPlot();


#Shallow lakes
Lk=c("ErieWest","Saginaw Bay","St. Clair"); 
#model comparison table
y=ModelComp(Lk,MaxKs=c(years=10,depths=10),type=2);
#Mussel biomass over depths and time
x=smoo2d2(k=20,nonzeroCenter=FALSE,rngMaxes=c(2018,16),plotgive=2,sampPt=list(pch=1,cex=0.5,col="black",jitter=0),dat=na.omit(cbind(log(1+ddall$DB),ddall$Year,ddall$Depth.m)[ddall$Lake%in%Lk,])); zCatPlot();
#Species composition
DomPlot(Lk);
#Mussel filtering intensity over depths and time
x=smoo2d2(k=10,nonzeroCenter=FALSE,rngMaxes=c(2018,16),n.grid=40,cats=1:3,plotgive=2,sampPt=c(pch=1,cex=0.5,col=1,jitter=0),dat=na.omit(cbind(Clear(),ddall$Year,ddall$Depth.m,ddall$Lake)[ddall$Lake%in%Lk,])); zCatPlot();




