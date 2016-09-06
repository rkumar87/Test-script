#
#
#


get_low_abundance_guides <- function(
	x,
	threshold=50,
	counts_col="total hits",
	genes_col='shrna.id'
	){
	return(x[which(x[,counts_col] < threshold),genes_col])
}


remove_controls <- function(x, pos="PLK1", neg="Olfr", id_colname="shrna.id"){
	ids <- x[,id_colname]
	pos_rows <-grep(pos,ids)
	neg_rows <-grep(neg,ids)
	pos_neg_rows <- c(pos_rows, neg_rows)
	return(x[-pos_neg_rows,])
}


get_pptm <- function(x){
	x_total_count_sum <- sum(
		x$total.hits
		)
	pptm <- x$total.hits / (x_total_count_sum / 10^7)
	pptm_psuedo <- make_pseudo_counts(pptm)
	x_with_pptm <- cbind(
		x,
		pptm_psuedo
		)
	return(x_with_pptm)
}

make_pseudo_counts <- function(x){
	pseudo_counts <- x + 0.5
	return(pseudo_counts)
}


# log2
# find ratio
# center and scale

znorm <- function(
	x0,
	x1,
	response_colname="pptm_psuedo",
	id_colname="shrna.id",
	min_counts=50,
	original_counts_colname="total.hits",
	exclude_guides
	){

	
	
	# find rows where the original counts for the
	# t0 group are below a threshold. exclude these
	# from both the t0 and t1 columns. Alternatively,
	# if a list of guides is provided (exclude_guides)
	# then use these to find rows_to_drop
	rows_to_drop <- NULL
	if(!is.null(exclude_guides)){
		rows_to_drop <- which(
			x0[,id_colname] %in% exclude_guides
			)
	}else{
		rows_to_drop <- which(
			x0[,original_counts_colname] < min_counts
			)
	}
	

	if(length(rows_to_drop) == 0){
		x0_cleaned <- x0
		x1_cleaned <- x1
	}else{
		x0_cleaned <- x0[-rows_to_drop,]
		x1_cleaned <- x1[-rows_to_drop,]
	}
	# log2 x0 and x1
	x0_log <- log2(x0_cleaned[,response_colname])
	x1_log <- log2(x1_cleaned[,response_colname])
	
	# find ratio
	x1_x0 <- x1_log / x0_log
	
	# center and scale
	x1_x0_med <- median(x1_x0)
	x1_x0_mad <- mad(x1_x0)
	x1_x0_zscore <- ( x1_x0 - x1_x0_med ) / x1_x0_mad
	names(x1_x0_zscore) <- x0_cleaned[,id_colname]
	tables_with_zscores <- cbind(
		x0_cleaned,
		x1_cleaned,
		x1_x0_zscore
		)
	
	dropped_data <- cbind(
		x0[rows_to_drop,],
		x1[rows_to_drop,],
		x1_x0_zscore=rep(NA, times=length(rows_to_drop))
		)
	tables_with_zscores <- rbind(
		tables_with_zscores,
		dropped_data
		)

	return(tables_with_zscores)
}

make_scatter_plot <- function(
	x0,
	x1,
	gene="PARP1",
	filename="plot.pdf",
	x0name="library",
	x1name="sample",
	response_colname="pptm_psuedo",
	main=""
	){
	
	pdf(filename, width=5, height=5)

	plot(
		x0[,response_colname],
		x1[,response_colname],
		log="xy",
		xlab=x0name,
		ylab=x1name,
		main=main
		)

	rows_to_mark <- grep(gene, x0$shrna.id)

	for(row in rows_to_mark){
		points(
			x0[row,response_colname],
			x1[row,response_colname],
			pch=19,
			col="red"
			)
	}

	dev.off()
}


correct_viability <- function(
	de_zscore_file,
	ve_zscore_file,
	plot_file="ve_correction.pdf"
	){
	
#	setwd("/Users/jfrankum/Documents/Crispr_screens/CRISPR_D05/")
#	ve_zscore_file="Sample10_9_filtered_viability_effect_Zscores_1W_270604.txt"
#	de_zscore_file="Sample12_10_filtered_VX_drug_effect_Zscores_1W_270604.txt"
	
	de <- read.table(
		file=de_zscore_file,
		header=TRUE,
		sep="\t",
		stringsAsFactors=FALSE
		)
	ve <- read.table(
		file=ve_zscore_file,
		header=TRUE,
		sep="\t",
		stringsAsFactors=FALSE
		)

	# ve and de have IDs in columns called 'shrna.id'
	# Z-scores and in the last column called 'x1_x0_zscore'
	
	# find IDs in common
	common_ids <- intersect(
		de$shrna.id,
		ve$shrna.id
		)
	
	de_ve_zscores <- cbind(
		de_z=de$x1_x0_zscore[which(de$shrna.id %in% common_ids)],
		ve_z=ve$x1_x0_zscore[which(ve$shrna.id %in% common_ids)]
		)
	rownames(de_ve_zscores) <- common_ids
	
	de_ve_lm <- lm(
		de_ve_zscores[,"de_z"] ~ de_ve_zscores[,"ve_z"]
		)
	
	de_corrected <- (de_ve_zscores[,"de_z"] - de_ve_lm$coefficients["(Intercept)"]) - (de_ve_lm$coefficients["de_ve_zscores[, \"ve_z\"]"] * de_ve_zscores[,"ve_z"])
	names(de_corrected) <- common_ids
	
#	y = m x + c
#	corrected de = de - c * (ve * m)
	
	
	
	pdf(
		plot_file,
		height=5,
		width=5
		)
		plot(
			de_ve_zscores[,"de_z"] ~ de_ve_zscores[,"ve_z"],
			xlab="viability effect",
			ylab="drug effect",
			main="before correction"
			)
		abline(
			de_ve_lm,
			col="red"
			)
		plot(
			de_corrected ~ de_ve_zscores[,"ve_z"],
			xlab="viability effect",
			ylab="corrected drug effect",
			main="after correction"
			)
		abline(
			lm(de_corrected ~ de_ve_zscores[,"ve_z"]),
			col="red"
			)			
	dev.off()
	
	de_with_ve_and_corrected <- cbind(
		de,
		ve=rep(NA, times=nrow(de)),
		corrected_de=rep(NA, times=nrow(de))
		)
	rows_to_replace <- which(
		de_with_ve_and_corrected$shrna.id %in% common_ids
		)
	
	de_with_ve_and_corrected[rows_to_replace,"ve"] <- ve$x1_x0_zscore[which(ve$shrna.id %in% common_ids)]
	de_with_ve_and_corrected[rows_to_replace,"corrected_de"] <- de_corrected
	
	return(de_with_ve_and_corrected)
}


combine_viability_corrected_data <- function(
	x,
	outfile="combined_viability_corrected_data.txt",
	wk
	){
	out_table <- NULL
	
	first_table <- x[["viability"]]
	
	if(colnames(first_table)[7] != "total.hits"){
		stop("Expected 7th column of viability table to be 'total.hits'")
	}
	if(colnames(first_table)[9] != "pptm_psuedo"){
		stop("Expected 9th column of viability table to be 'pptm_psuedo'")
	}
	if(colnames(first_table)[16] != "total.hits"){
		stop("Expected 16th column of viability table to be 'total.hits'")
	}
	if(colnames(first_table)[18] != "pptm_psuedo"){
		stop("Expected 18th column of viability table to be 'pptm_psuedo'")
	}
	
	
	out_table <-  cbind(
		first_table$shrna.id,
		first_table$shrna.sequence,
		first_table[,7],
		first_table[,9],
		first_table[,16],
		first_table[,18]
		)
	colnames(out_table) <- c(
		"guide.id",
		"guide.sequence",
		"dmso.t0.total.hits",
		"dmso.t0.pptm",
		paste("dmso_wk",wk,"_total.hits", sep=""),
		paste("dmso_wk",wk,"_pptm", sep="")
		)
	
	for(screen in names(x)){
		
		if(screen == "viability"){
			next
		}
		
		y <- x[[screen]]
		z <- cbind(
			y$total.hits.1,
			y$pptm_psuedo.1,
			y$x1_x0_zscore,
			y$ve,
			y$corrected_de
			)
		colnames(z) <- c(
			paste(screen,".total.hits", sep=""),
			paste(screen,".pptm", sep=""),
			paste(screen,".de_zscore", sep=""),
			paste(screen,".ve_zscore", sep=""),
			paste(screen,".de_corrected", sep="")
			)

		out_table <-  cbind(
			out_table,
			z
			)
	}
	
	write.table(
		out_table,
		file=outfile,
		col.names=TRUE,
		row.names=FALSE,
		sep="\t",
		quote=FALSE
		)
	
	return(out_table)
	
}




