###function to project points
mollweidePoints<-function(tcs.data){
points.theta <- tcs.data[, "h.theta"]
points.phi <- tcs.data[, "h.phi"]

n <- length(points.theta)

# Edges of the tetrahedron, adjusted
vert.theta <- c(-3.1415, 3.1415, -1.047198, 1.047198, -2.617994)
vert.phi <- c(-0.3398369, -0.3398369, -0.3398369, -0.3398369, 1.5707963)

# Edges of the figure
edge.theta <- c(-pi, -pi, pi, pi)
edge.phi <- c(-pi / 2, pi / 2, -pi / 2, pi / 2)

# adjust points

points.theta <- ifelse(points.theta >= -0.5235988,
  points.theta - (150 / 180 * pi),
  points.theta + (210 / 180 * pi)
)

# radians to degrees
coords.theta <- c(edge.theta, vert.theta, points.theta) * 180 / pi
coords.phi <- c(edge.phi, vert.phi, points.phi) * 180 / pi

# map projection coordinates
mp <- mapproj::mapproject(coords.theta, coords.phi, projection = "mollweide")

mp.p.theta <- mp$x[-c(seq_len(9))]
mp.p.phi <- mp$y[-c(seq_len(9))]
return(list(mp.p.phi,mp.p.theta))
}




###function to make shaded regions
mollweidePoly<-function(tcs.data,conc){
points.theta <- tcs.data[, "h.theta"]
points.phi <- tcs.data[, "h.phi"]

n <- length(points.theta)

# Edges of the tetrahedron, adjusted
vert.theta <- c(-3.1415, 3.1415, -1.047198, 1.047198, -2.617994)
vert.phi <- c(-0.3398369, -0.3398369, -0.3398369, -0.3398369, 1.5707963)

# Edges of the figure
edge.theta <- c(-pi, -pi, pi, pi)
edge.phi <- c(-pi / 2, pi / 2, -pi / 2, pi / 2)

# adjust points

points.theta <- ifelse(points.theta >= -0.5235988,
  points.theta - (150 / 180 * pi),
  points.theta + (210 / 180 * pi)
)

# radians to degrees
coords.theta <- c(edge.theta, vert.theta, points.theta) * 180 / pi
coords.phi <- c(edge.phi, vert.phi, points.phi) * 180 / pi

# map projection coordinates
mp <- mapproj::mapproject(coords.theta, coords.phi, projection = "mollweide")

mp.p.theta <- mp$x[-c(seq_len(9))]
mp.p.phi <- mp$y[-c(seq_len(9))]

####convert to concave hull
outline <-concaveman(cbind(mp.p.phi,mp.p.theta),concavity=conc)
return(outline)
}


###set up projection region
par(mar = c(0, 0, 0, 0))
plot(0, 0,
  axes = FALSE, xlab = "", ylab = "", type = "n", frame.plot = FALSE,
  xlim = c(-2, 2), ylim = c(-1, 1)
)
mapproj::map.grid(c(-180, 180, -90, 90),nx=9,ny=9, labels = FALSE, col = "grey")


myint=0.001             
myrep=0.5
	  
#####define boundaries:

#########################phagomixo###################################
##mix-rightSide
Mph_r<-seq(0.500,1,by=myint)
#Mphoto<-sample(rep(seq(0.501,1,by=myint),myrep))
Mphoto_r<-rev(seq(0.500,1,by=myint))
#Mproto<-sample(rep(seq(0.501,1,by=myint),myrep))
Mproto_r<-rep(0.500,(myrep/myint+1))
M_r<-cbind(rep("mix",(myrep/myint+1)),Mph_r,Mphoto_r,Mproto_r)

##mix-leftSide
Mph_l<-seq(0.500,1,by=myint)
#Mphoto<-sample(rep(seq(0.501,1,by=myint),myrep))
Mphoto_l<-rep(0.500,(myrep/myint+1))
#Mproto<-sample(rep(seq(0.501,1,by=myint),myrep))
Mproto_l<-rev(seq(0.500,1,by=myint))
M_l<-cbind(rep("mix",(myrep/myint+1)),Mph_l,Mphoto_l,Mproto_l)

###mix_bottom
Mph_b<-rep(0.500,(myrep/myint+1))
#Mphoto<-sample(rep(seq(0.501,1,by=myint),myrep))
Mphoto_b<-seq(0.500,1,by=myint)
#Mproto<-sample(rep(seq(0.501,1,by=myint),myrep))
Mproto_b<-rev(seq(0.500,1,by=myint))
M_b<-cbind(rep("mix",(myrep/myint+1)),Mph_b,Mphoto_b,Mproto_b)


##############photoauto#############################################
###photoauto_rightSide
Aph_r<-seq(0,0.499,by=myint)
Aphoto_r<-seq(0.501,1,by=myint)
Aproto_r<-rep(0.4999,(myrep/myint))
A_r<-cbind(rep("photoauto",(myrep/myint)),Aph_r,Aphoto_r,Aproto_r)

###photoauto_leftSide
Aph_l<-seq(0,0.499,by=myint)
Aphoto_l<-rep(0.500,(myrep/myint))
Aproto_l<-seq(0.501,1,by=myint)
A_l<-cbind(rep("photoauto",(myrep/myint)),Aph_l,Aphoto_l,Aproto_l)

###photoauto_top
Aph_t<-rep(0.499,myrep/myint)
Aphoto_t<-seq(0.501,1,by=myint)
Aproto_t<-rev(seq(0.501,1,by=myint))
A_t<-cbind(rep("photoauto",myrep/myint),Aph_t,Aphoto_t,Aproto_t)

myint=0.0001
###photoauto_bottom
Aphoto_br<-seq(0.5,0.601,by=myint)
Aproto_br<-Aphoto_br-myint 
Aproto_bl<-seq(0.5,0.601,by=myint)
Aphoto_bl<-Aproto_bl-myint
Aph_br<-rep(0,length(c(Aphoto_br,Aphoto_bl)))
Aph_bl<-rep(0,length(c(Aphoto_br,Aphoto_bl)))
A_b<-cbind(rep("photoauto",length(c(Aphoto_br,Aphoto_bl))),c(Aph_br,Aph_bl),c(Aphoto_br,Aphoto_bl),c(Aproto_br,Aproto_bl))
myint=0.001


###################phago_proto########################################
##PHPR-leftSide
PHPROTOph_l<-seq(0.501,1,by=myint)
PHPROTOphoto_l<-seq(0,0.499,by=myint)
PHPROTOproto_l<-rep(0.501,myrep/myint)
PHPROTO_l<-cbind(rep("phagoosmo",myrep/myint),PHPROTOph_l,PHPROTOphoto_l,PHPROTOproto_l)

##PHPR-bottom
PHPROTOph_b<-rep(0.501,myrep/myint)
PHPROTOphoto_b<-seq(0,0.499,by=myint)
PHPROTOproto_b<-seq(0.501,1,by=myint)
PHPROTO_b<-cbind(rep("phagoosmo",myrep/myint),PHPROTOph_b,PHPROTOphoto_b,PHPROTOproto_b)

##PHPR-rightSide
PHPROTOph_r<-seq(0.501,1,by=myint)
PHPROTOphoto_r<-rep(0.499,myrep/myint)
PHPROTOproto_r<-rev(seq(0.501,1,by=myint))
PHPROTO_r<-cbind(rep("phagoosmo",myrep/myint),PHPROTOph_r,PHPROTOphoto_r,PHPROTOproto_r)



###################photo_only########################################
##PHOTOonly-leftSide
PHOTOonlyph_l<-seq(0,0.499,by=myint)
PHOTOonlyphoto_l<-seq(0.501,1,by=myint)
PHOTOonlyproto_l<-rep(0.4989,myrep/myint)
PHOTOonly_l<-cbind(rep("photoonly",myrep/myint),PHOTOonlyph_l,PHOTOonlyphoto_l,PHOTOonlyproto_l)

##PHOTOonly-top
PHOTOonlyph_t<-rep(0.499,myrep/myint)
PHOTOonlyphoto_t<-seq(0.501,1,by=myint)
PHOTOonlyproto_t<-seq(0,0.499,by=myint)
PHOTOonly_t<-cbind(rep("photoonly",myrep/myint),PHOTOonlyph_t,PHOTOonlyphoto_t,PHOTOonlyproto_t)

##PHOTOonly-rightSide
PHOTOonlyph_r<-rev(seq(0,0.499,by=myint))
PHOTOonlyphoto_r<-rep(0.501,myrep/myint)
PHOTOonlyproto_r<-seq(0,0.499,by=myint)
PHOTOonly_r<-cbind(rep("photoonly",myrep/myint),PHOTOonlyph_r,PHOTOonlyphoto_r,PHOTOonlyproto_r)


###################proto_only########################################
##PROTOonly-rightSide
PROTOonlyph_r<-seq(0,0.499,by=myint)
PROTOonlyphoto_r<-rep(0.499,myrep/myint)
PROTOonlyproto_r<-seq(0.501,1,by=myint)
PROTOonly_r<-cbind(rep("protoonly",myrep/myint),PROTOonlyph_r,PROTOonlyphoto_r,PROTOonlyproto_r)

##PROTOonly-top
PROTOonlyph_t<-rep(0.4989,myrep/myint)
PROTOonlyphoto_t<-seq(0,0.499,by=myint)
PROTOonlyproto_t<-seq(0.501,1,by=myint)
PROTOonly_t<-cbind(rep("protoonly",myrep/myint),PROTOonlyph_t,PROTOonlyphoto_t,PROTOonlyproto_t)

##PROTOonly-leftSide
PROTOonlyph_l<-seq(0,0.499,by=myint)
PROTOonlyphoto_l<-rev(seq(0,0.499,by=myint))
PROTOonlyproto_l<-rep(0.501,myrep/myint)
PROTOonly_l<-cbind(rep("protoonly",myrep/myint),PROTOonlyph_l,PROTOonlyphoto_l,PROTOonlyproto_l)


###################phago_only########################################
##PHAGOonly-rightSide
myint=0.0001
PHAGOonlyph_r<-seq(0.5001,1,by=myint)
PHAGOonlyphoto_r<-rep(0.4999,myrep/myint)
PHAGOonlyproto_r<-seq(0,0.4999,by=myint)
PHAGOonly_r<-cbind(rep("phagoonly",myrep/myint),PHAGOonlyph_r,PHAGOonlyphoto_r,PHAGOonlyproto_r)

myint=0.001
##PHAGOonly-bottom
PHAGOonlyph_b<-rep(0.501,myrep/myint)
PHAGOonlyphoto_b<-seq(0,0.499,by=myint)
PHAGOonlyproto_b<-rev(seq(0,0.499,by=myint))
PHAGOonly_b<-cbind(rep("phagoonly",myrep/myint),PHAGOonlyph_b,PHAGOonlyphoto_b,PHAGOonlyproto_b)

##PHAGOonly-leftSide
PHAGOonlyph_l<-seq(0.501,1,by=myint)
PHAGOonlyphoto_l<-seq(0,0.499,by=myint)
PHAGOonlyproto_l<-rep(0.499,myrep/myint)
PHAGOonly_l<-cbind(rep("phagoonly",myrep/myint),PHAGOonlyph_l,PHAGOonlyphoto_l,PHAGOonlyproto_l)

##PHAGOonly-top
PHAGOonlyphoto_tr<-seq(0,0.401,by=myint)
PHAGOonlyproto_tr<-PHAGOonlyphoto_tr-myint
PHAGOonlyproto_tl<-seq(0,0.401,by=myint)
PHAGOonlyphoto_tl<-PHAGOonlyproto_tl-myint
PHAGOonlyph_tr<-rep(1,length(PHAGOonlyphoto_tr))
PHAGOonlyph_tl<-rep(1,length(PHAGOonlyphoto_tl))
PHAGOonly_tt<-cbind(rep("phagoonly",length(c(PHAGOonlyphoto_tr,PHAGOonlyphoto_tl))),c(PHAGOonlyph_tr,PHAGOonlyph_tl),c(PHAGOonlyphoto_tl,PHAGOonlyphoto_tr),c(PHAGOonlyproto_tl,PHAGOonlyproto_tr))

PHAGOonlyphoto_tr<-seq(0,0.401,by=myint)
PHAGOonlyproto_tr<-PHAGOonlyphoto_tr-myint
PHAGOonlyproto_tl<-seq(0,0.401,by=myint)
PHAGOonlyphoto_tl<-PHAGOonlyproto_tl-myint
PHAGOonlyph_tr<-sample(seq(0.501,1,by=myint),length(PHAGOonlyphoto_tr))
PHAGOonlyph_tl<-sample(seq(0.501,1,by=myint),length(PHAGOonlyphoto_tl))
PHAGOonly_tb<-cbind(rep("phagoonly",length(c(PHAGOonlyphoto_tr,PHAGOonlyphoto_tl))),c(PHAGOonlyph_tr,PHAGOonlyph_tl),c(PHAGOonlyphoto_tl,PHAGOonlyphoto_tr),c(PHAGOonlyproto_tl,PHAGOonlyproto_tr))

###################phagophoto########################################

##PHPHOTO-leftSide
PHPHOTOph_l<-seq(0.501,1,by=myint)
PHPHOTOphoto_l<-rep(0.501,myrep/myint)
PHPHOTOproto_l<-seq(0,0.499,by=myint)
PHPHOTO_l<-cbind(rep("phagophoto",myrep/myint),PHPHOTOph_l,PHPHOTOphoto_l,PHPHOTOproto_l)

##PHPHOTO-bottom
PHPHOTOph_b<-rep(0.501,myrep/myint)
PHPHOTOphoto_b<-seq(0.501,1,by=myint)
PHPHOTOproto_b<-seq(0,0.499,by=myint)
PHPHOTO_b<-cbind(rep("phagophoto",myrep/myint),PHPHOTOph_b,PHPHOTOphoto_b,PHPHOTOproto_b)

##PHPHOTO-rightSide
PHPHOTOph_r<-seq(0.501,1,by=myint)
PHPHOTOphoto_r<-rev(seq(0.501,1,by=myint))
PHPHOTOproto_r<-rep(0.499,myrep/myint)
PHPHOTO_r<-cbind(rep("phagophoto",myrep/myint),PHPHOTOph_r,PHPHOTOphoto_r,PHPHOTOproto_r)


###################parasite########################################

##PARA-rightleftSide
PARAph_rl<-seq(0.001,0.499,by=myint)
PARAphoto_rl<-rep(0.499,(myrep/myint-1))
PARAproto_rl<-rev(seq(0.001,0.499,by=myint))
PARA_rl<-cbind(rep("para_r",(myrep/myint-1)),PARAph_rl,PARAphoto_rl,PARAproto_rl)

##PARA-top
PARAph_t<-rep(0.499,myrep/myint)
PARAphoto_t<-seq(0,0.499,by=myint)
PARAproto_t<-rev(seq(0,0.499,by=myint))
PARA_t<-cbind(rep("para",myrep/myint),PARAph_t,PARAphoto_t,PARAproto_t)
PARA_tl<-PARA_t[1:250,]
PARA_tl[,1]<-rep("para_l",length(PARA_tl[,1]))
PARA_tr<-PARA_t[251:500,]
PARA_tr[,1]<-rep("para_r",length(PARA_tr[,1]))

##PARA-leftrightSide
PARAph_lr<-seq(0,0.499,by=myint)
PARAphoto_lr<-rev(seq(0,0.499,by=myint))
PARAproto_lr<-rep(0.499,myrep/myint)
PARA_lr<-cbind(rep("para_l",myrep/myint),PARAph_lr,PARAphoto_lr,PARAproto_lr)

##PARA-edges
PARAph_ll1<-seq(0,0.499,by=myint)
PARAphoto_ll1<-rep(0.249,myrep/myint)
PARAproto_ll1<-rep(0.25,myrep/myint)
PARA_ll1<-cbind(rep("para_l",myrep/myint),PARAph_ll1,PARAphoto_ll1,PARAproto_ll1)

myints<-0.00001
PARAph_ll2<-seq(0,0.499,by=myints)
PARAproto_ll2<-rep(0.499,(myrep/myints-98))
PARAphoto_ll2<-PARAproto_ll2-myints
PARA_ll2<-cbind(rep("para_l",length(PARAproto_ll2)),PARAph_ll2,PARAphoto_ll2,PARAproto_ll2)

PARAph_rr1<-seq(0,0.499,by=myint)
PARAphoto_rr1<-rep(0.25,myrep/myint)
PARAproto_rr1<-rep(0.249,myrep/myint)
PARA_rr1<-cbind(rep("para_r",myrep/myint),PARAph_rr1,PARAphoto_rr1,PARAproto_rr1)

myints<-0.00001
PARAph_rr2<-seq(0,0.499,by=myints)
#PARAph_rr2<-PARAph_rr2[1:(i+100)==(i+100)]
PARAphoto_rr2<-rep(0.499,(myrep/myints-98))
#PARAphoto_rr2<-PARAphoto_rr2[1:(i+100)==(i+100)]
PARAproto_rr2<-PARAphoto_rr2-myints
PARA_rr2<-cbind(rep("para_r",length(PARAproto_rr2)),PARAph_rr2,PARAphoto_rr2,PARAproto_rr2)



##########################################################################################
#combine boundaries
##########################################################################################

std.df<-rbind(M_r,M_l,M_b,A_r,A_l,A_t,A_b,PHPROTO_r,PHPROTO_l,PHPROTO_b,PHOTOonly_l,PHOTOonly_t,PHOTOonly_r,PROTOonly_l,PROTOonly_t,PROTOonly_r,PHAGOonly_l,PHAGOonly_r,PHAGOonly_b,PHAGOonly_tt,PHAGOonly_tb,PHPHOTO_l,PHPHOTO_r,PHPHOTO_b,PARA_rl,PARA_tl,PARA_tr,PARA_lr,PARA_ll1,PARA_ll2,PARA_rr1,PARA_rr2)
std.df<-as.data.frame(std.df)
std.df[,2:4] = map_df(std.df[,2:4], as.character)
std.df[,2:4] = map_df(std.df[,2:4], as.numeric)
std.df$nopred <- 1-(rowSums(std.df[2:4]))/3

col.df<-std.df[std.df[,1]=="mix",2:5]
colnames(col.df)<-c("u","s","m","l")
tcs.mix<-tcspace(col.df)
mix.pts<-mollweidePoly(tcs.mix,1)

col.df<-std.df[std.df[,1]=="photoauto",2:5]
colnames(col.df)<-c("u","s","m","l")
tcs.mix<-tcspace(col.df)
photoauto.pts<-mollweidePoly(tcs.mix,1)


col.df<-std.df[std.df[,1]=="phagoosmo",2:5]
colnames(col.df)<-c("u","s","m","l")
tcs.mix<-tcspace(col.df)
phagoosmo.pts<-mollweidePoly(tcs.mix,1)

col.df<-std.df[std.df[,1]=="photoonly",2:5]
colnames(col.df)<-c("u","s","m","l")
tcs.mix<-tcspace(col.df)
photoonly.pts<-mollweidePoly(tcs.mix,1)

col.df<-std.df[std.df[,1]=="protoonly",2:5]
colnames(col.df)<-c("u","s","m","l")
tcs.mix<-tcspace(col.df)
osmo.pts<-mollweidePoly(tcs.mix,1)

col.df<-std.df[std.df[,1]=="phagoonly",2:5]
colnames(col.df)<-c("u","s","m","l")
tcs.mix<-tcspace(col.df)
phago.pts<-mollweidePoly(tcs.mix,1)


col.df<-std.df[std.df[,1]=="phagophoto",2:5]
colnames(col.df)<-c("u","s","m","l")
tcs.mix<-tcspace(col.df)
phagophoto.pts<-mollweidePoly(tcs.mix,1)

col.df<-std.df[std.df[,1]=="para_l",2:5]
colnames(col.df)<-c("u","s","m","l")
tcs.mix<-tcspace(col.df)
para_l.pts<-mollweidePoly(tcs.mix,1)

col.df<-std.df[std.df[,1]=="para_r",2:5]
colnames(col.df)<-c("u","s","m","l")
tcs.mix<-tcspace(col.df)
para_r.pts<-mollweidePoly(tcs.mix,1)

####plot polygons of shaded regions
polygon(mix.pts[,1] ~ mix.pts[,2],col=rgb(0,0,1,0.25),lty=0)
polygon(photoauto.pts[,1] ~ photoauto.pts[,2],col=rgb(0,1,0,0.25),lty=0)
polygon(phagoosmo.pts[,1] ~ phagoosmo.pts[,2],col=rgb(0.7254902,0.4980392,0.4980392,0.25),lty=0)
polygon(photoonly.pts[,1] ~ photoonly.pts[,2],col=rgb(0.7254902,0.9333333,0.6745098,0.25),lty=0)
polygon(osmo.pts[,1] ~ osmo.pts[,2],col=rgb(0.7529412,0.7529412,0.7529412,0.25),lty=0)
polygon(phago.pts[,1] ~ phago.pts[,2],col=rgb(1,0,0,0.25),lty=0)
polygon(phagophoto.pts[,1] ~ phagophoto.pts[,2],col=rgb(1,0.8235294,0.5843137,0.25),lty=0)
polygon(para_l.pts[,1] ~ para_l.pts[,2],col=rgb(0.3529412,0.3529412,0.3529412,0.25),lty=0)
polygon(para_r.pts[,1] ~ para_r.pts[,2],col=rgb(0.3529412,0.3529412,0.3529412,0.25),lty=0)


