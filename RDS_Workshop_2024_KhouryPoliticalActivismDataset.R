##################################################################
## R code - used to analyze the Khoury Syrian political activism dataset
##
## Author: Chris Meaney
## Date: January 2024
###########################################################

####################
## Package dependencies
####################
library(RDS)


####################
## Import the Khoury political activism dataset into the R environment
####################
fpath <- "C://Users//ChristopherMeaney//Desktop//tmp//RDS_Workshop//Winter2024//ZZZ_RDS_ExampleDataset//dataverse_files//rds_replication_data.csv"
X <- read.csv(file=fpath, header=TRUE, sep=",", stringsAsFactors=FALSE)

## Some properties of the imported dataFrame
data.frame(names_=names(X), class_=sapply(X, class))

## Uniqueness of IDs
nrow(X)
length(unique(X$id))
length(unique(X$recruit.id))

## Drop variables that will not be used in this analysis
keep_vars <- c("id",
				"recruit.id",
				"coupon.1",
				"coupon.2",
				"coupon.3",
				"degree",
				"sex",
				"age",
				"syria_pre2011",
				"syria_post2011")

X <- X[, names(X) %in% keep_vars]
dim(X)


#################################
## Recodes of demographic and outcome variables
#################################

## Age
summary(X$age)
X$age <- as.numeric(X$age)
## Sex 
table(X$sex)
X$sex_R <- as.factor(ifelse(X$sex==1, "Male", "Female"))
table(X$sex_R, X$sex)
data.frame(table(X$sex_R))

## Involvement in political activism pre-2011
table(X$syria_pre2011)
X$syria_pre <- as.factor(ifelse(X$syria_pre2011==1, "Yes", "No"))
table(X$syria_pre, X$syria_pre2011, useNA="always")
data.frame(table(X$syria_pre, useNA="always"))

## Involvement in political activism post-2011
table(X$syria_post2011)
X$syria_post <- as.factor(ifelse(X$syria_post==1, "Yes", "No"))
table(X$syria_post, X$syria_post2011, useNA="always")
data.frame(table(X$syria_post, useNA="always"))


###################
## Inspect some properties of the RDS design, create an RDS dataFrame object, and use this object to generate some RDS diagnostic plots (and eventually estimators)
###################

## Use built-in RDS functionality to determine recruiter.id from various ID and coupon info
X$recruiter.id <- rid.from.coupons(data=X,
									 subject.coupon="recruit.id",
									 coupon.variables=paste0("coupon.", 1:3),
									 subject.id="id")

## Flag for seeds
X$seed <- ifelse(X$recruiter.id=="seed", "seed", "non-seed")
table(X$seed)

## Flag for wave
X$wave <- as.numeric(nchar(X$recruit.id)) - 1
table(X$wave)

## Pop size estimates
pop_size <- 50000

## Create an RDS dataFrame --- as special data structure needed to enable certain functionality in RDS package
X_rds <- as.rds.data.frame(df=X, 
							id="id", 
							recruiter.id="recruiter.id", 
							network.size="degree", 
							time="wave",
							population.size=pop_size, 
							max.coupons=3,
							notes=NULL
							)

## Assert whether the created RDS dataFrame is valid
assert.valid.rds.data.frame(X_rds)

## Plot a recruitment tree for persons enrolled in this RDS study
plot(x=X_rds, 
	plot.type="Recruitment tree", 
	stratify.by="seed",
	main="RDS Recruitment Tree")

##
## ?plot.rds.data.frame
##
## Different types of diagnostic plots which can be created by calling the "plot" function against an RDS dataFrame
##
## "Recruitment tree"
## "Network size by wave"
## "Recruits by wave",
## "Recruits per seed"
## "Recruits per subject"
##

##
## Some other RDS diagnostic plots...
##
plot(x=X_rds, plot.type="Network size by wave", main="Network Size by Wave")
plot(x=X_rds, plot.type="Recruits by wave", main="Recruits by Wave")
plot(x=X_rds, plot.type="Recruits per seed", main="Recruits per seed")
plot(x=X_rds, plot.type="Recruits per subject", main="Recruits per subject")



################################################
## Create two new function to generate comparative RDS estimates of categorical/numeric variables 
################################################

## Function to compare estimated proportions from traditional methods vs RDS methods
get_rds_binary_props <- function(dat, varx, n) {

	## Simple counts and percentages
	var_names <- names(dat)
	xdf <- dat[!(is.na(dat[[varx]])), ]

	counts <- table(xdf[[varx]])
	percents <- prop.table(counts)
	counts_percents <- cbind(counts, 
							#(percents*100, 2)
							percents
							)

	## Exact binomial estimator
	binom_test <- binom.test(counts[2], sum(counts))
	binom_est <- binom_test$estimate
	binom_ll <- binom_test$conf.int[1]
	binom_ul <- binom_test$conf.int[2]
	binom_width <- binom_ul - binom_ll
	binom_vec <- c(binom_est, binom_ll, binom_ul, binom_width)
	binom_df <- data.frame(ExactBinomial=binom_vec)

	## RDS-1 estimator
	rds1est <- RDS.I.estimates(rds.data=xdf, outcome.variable=varx)
	rds1ints <- data.frame(rds1est$interval)
	rds1ints$ci_width <- rds1ints$upper - rds1ints$lower
	rds1ints_sm <- data.frame(t(rds1ints[1:length(counts), c(1:3,7)]))
	names(rds1ints_sm) <- paste0(rownames(rds1ints), "_RDS1")
	# rds1ints_sm

	## RDS-2 estimator
	rds2ests <- RDS.II.estimates(rds.data=xdf, outcome.variable=varx)
	rds2ints <- data.frame(rds2ests$interval)
	rds2ints$ci_width <- rds2ints$upper - rds2ints$lower
	rds2ints_sm <- data.frame(t(rds2ints[1:length(counts), c(1:3,7)]))
	names(rds2ints_sm) <- paste0(rownames(rds2ints), "_RDS2")
	# rds2ints_sm

	## Gile's estimator
	rdsSSests <- RDS.SS.estimates(rds.data=xdf, outcome.variable=varx, N=n)
	rdsSSints <- data.frame(rdsSSests$interval)
	rdsSSints$ci_width <- rdsSSints$upper - rdsSSints$lower
	rdsSSints_sm <- data.frame(t(rdsSSints[1:length(counts), c(1:3,7)]))
	names(rdsSSints_sm) <- paste0(rownames(rdsSSints), "_RDS_SS")
	# rdsSSints_sm

	## HCG estimator
	rdsHCGests <- RDS.HCG.estimates(rds.data=xdf, outcome.variable=varx, N=n)
	rdsHCGints <- data.frame(rdsHCGests$interval)
	rdsHCGints$ci_width <- rdsHCGints$upper - rdsHCGints$lower
	rdsHCGints_sm <- data.frame(t(rdsHCGints[1:length(counts), c(1:3,7)]))
	names(rdsHCGints_sm) <- paste0(rownames(rdsHCGints), "_RDS_HCG")
	# rdsHCGints_sm

	## Bind the RDS estimators togethers
	rds_bind_df <- data.frame(cbind(binom_df,
								rds1ints_sm[,2], 
								rds2ints_sm[,2], 
								rdsSSints_sm[,2],
								rdsHCGints_sm[,2]))
	
	names(rds_bind_df) <- c("ExactBinomial","RDS1","RDS2","RDS_SS","RDS_HCG")
	rownames(rds_bind_df) <- c("PointEst", "LL95CI", "UL95CI","CI_Width")

	## Return the dataFrame to the user
	return(rds_bind_df)

}


## Function to compare estimated means using traditional methods versus RDS methods
get_rds_numeric <- function(dat, varx, n) {

	## Simple counts and percentages
	xdf <- dat[!(is.na(dat[[varx]])), ]

	n <- nrow(xdf)
	mu <- mean(xdf[[varx]], na.rm=TRUE)
	se <- sqrt(var(xdf[[varx]], na.rm=TRUE)/n)

	## T-test
	t_test <- t.test(xdf[[varx]])
	t_est <- t_test$estimate
	t_ll <- t_test$conf.int[1]
	t_ul <- t_test$conf.int[2]
	t_width <- t_ul - t_ll
	t_vec <- c(t_est, t_ll, t_ul, t_width)
	t_df <- data.frame(t_test=t_vec)

		## RDS-1 estimator
	rds1est <- RDS.I.estimates(rds.data=xdf, outcome.variable=varx)
	rds1ints <- data.frame(rds1est$interval)
	rds1ints$ci_width <- rds1ints$upper - rds1ints$lower
	rds1ints_sm <- data.frame(t(rds1ints[1, c(1:3,7)]))
	names(rds1ints_sm) <- paste0(rownames(rds1ints), "_RDS1")
	# rds1ints_sm

	## RDS-2 estimator
	rds2ests <- RDS.II.estimates(rds.data=xdf, outcome.variable=varx)
	rds2ints <- data.frame(rds2ests$interval)
	rds2ints$ci_width <- rds2ints$upper - rds2ints$lower
	rds2ints_sm <- data.frame(t(rds2ints[1, c(1:3,7)]))
	names(rds2ints_sm) <- paste0(rownames(rds2ints), "_RDS2")
	# rds2ints_sm

	## Gile's estimator
	rdsSSests <- RDS.SS.estimates(rds.data=xdf, outcome.variable=varx, N=n)
	rdsSSints <- data.frame(rdsSSests$interval)
	rdsSSints$ci_width <- rdsSSints$upper - rdsSSints$lower
	rdsSSints_sm <- data.frame(t(rdsSSints[1, c(1:3,7)]))
	names(rdsSSints_sm) <- paste0(rownames(rdsSSints), "_RDS_SS")
	# rdsSSints_sm

	## HCG estimator
	rdsHCGests <- RDS.HCG.estimates(rds.data=xdf, outcome.variable=varx, N=n)
	rdsHCGints <- data.frame(rdsHCGests$interval)
	rdsHCGints$ci_width <- rdsHCGints$upper - rdsHCGints$lower
	rdsHCGints_sm <- data.frame(t(rdsHCGints[1, c(1:3,7)]))
	names(rdsHCGints_sm) <- paste0(rownames(rdsHCGints), "_RDS_HCG")
	# rdsHCGints_sm

	## Bind the RDS estimators togethers
	rds_bind_df <- data.frame(cbind(t_df,
								rds1ints_sm, 
								rds2ints_sm, 
								rdsSSints_sm,
								rdsHCGints_sm))
	
	names(rds_bind_df) <- c("T_test","RDS1","RDS2","RDS_SS","RDS_HCG")
	rownames(rds_bind_df) <- c("PointEst", "LL95CI", "UL95CI","CI_Width")

	## Return the dataFrame to the user
	return(rds_bind_df)

}


################################################
## Compute and compare unadjusted vs. adjusted estimates of outcome variables (political activism) using RDS weighting methods
################################################


##
## Involvement in political activism pre-2011
##
get_rds_binary_props(dat=X_rds, varx="syria_pre", n=50000)

##
## Involvement in political activism post-2011
##
get_rds_binary_props(dat=X_rds, varx="syria_post", n=50000)

##
## Sex (binary)
##
get_rds_binary_props(dat=X_rds, varx="sex_R", n=50000)

##
## Age (continuous)
##
get_rds_numeric(dat=X_rds, varx="age", n=50000)



#################################################
##
## Inference --- does distribution of post-syrian activism differ as a function of demographic variables (e.g. sex)
## 
## A newish function bootstrap.contingency.test() allows for such an inference 
## Estimates HCG pop weighted estimated of chi-square test between two categorical variables
##
#################################################

##
## Are men and women equally likely to participate in political activism post-2011; or are men != women
##

## Contingency table
table(x=X_rds$syria_post, y=X_rds$sex_R)
## Straitfied proportions
prop.table(table(x=X_rds$syria_post, y=X_rds$sex_R), margin=2)
## Chi-square test
chisq.test(table(X_rds$syria_post, X_rds$sex_R))
chisq.test(table(X_rds$sex_R, X_rds$syria_post))
## Or Wald Z-test on independent binomial proportions
prop.test(table(X_rds$syria_post, X_rds$sex_R))
prop.test(table(X_rds$sex_R, X_rds$syria_post))

## Bootstrap HCG chi-square test from the RDS package
bootstrap.contingency.test(rds.data=X_rds, row.var="syria_post", col.var="sex_R", weight.type="HCG", number.of.bootstrap.samples=50, verbose=TRUE)
bootstrap.contingency.test(rds.data=X_rds, row.var="sex_R", col.var="syria_post", weight.type="HCG", number.of.bootstrap.samples=50, verbose=TRUE)

## NOTE: the lack of symmetry of the test above

##
## Other weighting options
##
## "HCG"
## "RDS-II"
## "Arithmetic Mean"
##


##############################################
##
## Bottleneck plots (for assessing convergence of RDS estimates)
##
##############################################


## Bottleneck plot illustrating convergence of estimate regarding the proportion of persons participating in syrian activism pre-2011 (across different seeds)
prop.table(table(X_rds$syria_pre))
bottleneck.plot(X_rds, "syria_pre")

## Bottleneck plot illustrating convergence of estimate regarding the proportion of persons participating in syrian activism post-2011 (across different seeds)
prop.table(table(X_rds$syria_post))
bottleneck.plot(X_rds, "syria_post")


##############################################
##
## RDS convergence plots
##
##############################################
convergence.plot(rds.data=X_rds,
					outcome.variable="syria_post",
					est.func=RDS.II.estimates,
					as.factor=FALSE,
					n.eval.points=25
					)


##############################################
##
## Compute RDS weights (based on various weighting schemes), which can then be used in different RDS estimators --- for example, RDS regression
##
##############################################

hcg_weights_syria_post <- compute.weights(rds.data=X_rds,
										outcome.variable="syria_post",
										weight.type="HCG",
										N=50000,
										subset=NULL,
										control=control.rds.estimates()
										)

##
## RDS weight options 
##
## c("Gile's SS", "RDS-I", "RDS-I (DS)", "RDS-II", "Arithmetic Mean", "HCG")
##

##
## Compare unweighted vs weighted RDS regression
##
syria_post_reg <- glm(I(syria_post=="Yes") ~ sex_R + age, data=X_rds, family="binomial"(link="logit"))
summary(syria_post_reg)

syria_post_reg_weighted <- glm(I(syria_post=="Yes") ~ sex_R + age, data=X_rds, family="binomial"(link="logit"), weights=hcg_weights_syria_post)
summary(syria_post_reg_weighted)


##############################################
##
## Differential activity estimates --- compares "mean degree estimates" in one group vs another group
##
##############################################
differential.activity.estimates(rds.data=X_rds,
								outcome.variable="syria_post",
								weight.type = "HCG",
								N = 50000,
								subset = NULL
								)


################################################
##
##  get.net.size() returns degree of members in RDS design
##
################################################
get.net.size(X_rds)
summary(get.net.size(X_rds))


################################################
##
##  get.number.of.recruits() returns number of recruits for each member of the RDS design
##
################################################
get.number.of.recruits(X_rds)
table(get.number.of.recruits(X_rds))


#################################################
##
## Homophily estimates 
##
#################################################
homophily.estimates(rds.data=X_rds,
					outcome.variable="syria_post",
					weight.type = "HCG",
					uncertainty = NULL,
					recruitment = FALSE,
					N = 50000,
					to.group0.variable = NULL,
					to.group1.variable = NULL,
					number.ss.samples.per.iteration = NULL,
					confidence.level = 0.95
					)



