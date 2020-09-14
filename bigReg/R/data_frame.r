#' @name plasma
#' @title plasma data from the HSAUR package
#' @description Dataset from the HSAUR package
#' @usage data(plasma)
#' @docType data
#'
#' @format a data.frame
#'
#' @keywords datasets
#'
#' @references HSAUR R package
#' (\href{https://cran.r-project.org/package=HSAUR}{HSAUR package})
#' @source \href{https://cran.r-project.org/package=HSAUR}{HSAUR package}
#' @details ...
#'
#' @examples
#' data(plasma)
#' head(plasma)
#' @export
NULL

#' Creates a data_frame object
#' @title function to create a data_frame object
#' @description function to create a data_frame object. The data_frame object
#'              is an object that is held on disk. It is written to a folder
#'              \code{path} on disk where the data is written to in blocks or
#'              chunks. The data is written in binary format using a C++ function
#'              in purely numerical data and a mapping to the table is held in
#'              a ".meta_data" file in the folder. The table object accomodates
#'              numeric, factor, and character (converted to factor).
#' @param data data.frame object to be converted into a data_frame object
#' @param chunkSize number of rows to be used in each chunk
#' @param path character to folder where the object will be created
#' @param nCores the number of cores to use defaults to parallel::detectCores()
#' @param ... not currently used.
#' @examples
#' irisA <- data_frame(iris[1:75,], 10, "irisA", nCores = 1)
#' irisA$append(iris[76:150,])
#' irisA$head()
#' irisA$tail(10)
#' irisA$delete(); rm(irisA)
#' @export
data_frame <- function(data = stop("data must be supplied"),
	chunkSize = stop("chunkSize must be specified, a good number is 50000"),
	path = stop("path must be specified"), nCores = parallel::detectCores(), ...){

	this <- new.env()
	this$size <- vector(mode = "numeric", length = 0)

	if(!dir.exists(path)){
		dir.create(path)
	}

	this$path <- path

	# Default lapply function
	# this$lApply <- lapply
	this$lApply <- function(x, fun, ...)mclapply(x, fun, ..., mc.cores = nCores)

	set_lapply_function <- function(fun){
		this$lApply <- fun
	}
	environment(set_lapply_function) <- this
	this$set_lapply_function <- set_lapply_function
	
	this$chunkSize <- chunkSize
	this$path <- path

	this$isFile <- function(path){
		is(path, "character") && (!dir.exists(path) & file.exists(path))
	}

	#this$isFile <- isFile

	set_reader <- function(fun){
		environment(fun) <- this
		this$Reader <- fun
	}
	environment(set_reader) <- this
	this$set_reader <- set_reader

	this$fReader <- read.csv

	set_file_reader <- function(fun){
		environment(fun) <- this
		this$fReader <- fun
	}

	environment(set_file_reader) <- this
	this$set_file_reader <- set_file_reader

	# This reads the file if it is a data path
	if(this$isFile(data)){
		data <- this$fReader(data)
	}

	get_unique_values <- function(col_name = stop("col_name: has not been specified"), 
										data = stop("data: has not been specified")){
		column <- data[, col_name]
		if(is(column[1], "factor")){
			uniques <- levels(column[1])
		}else if(is(column[1], "character")){
			uniques <- unique(column)
		}else{
			stop("The data submitted is neither factor or character")
		}
		return(uniques)
	}

	this$get_unique_values <- get_unique_values


	# Function to read colNames from data frame
	get_df_col_classes <- function(data = stop("No Data Submitted!")){
		vapply(data[1,], class, FUN.VALUE = "")
	}

	this$get_df_col_classes <- get_df_col_classes

	this$col_names <- names(data)
	col_classes <- this$get_df_col_classes(data)

	this$ncols <- ncol(data)
	this$nrows <- 0
	this$list_frame <- vector(mode = "list", length = this$ncols)
	names(this$list_frame) <- this$col_names
	factor_cols <- names(col_classes[which(col_classes %in% c("character", "factor"))])

	this$factors <- list()
	if(length(factor_cols) > 0){
		this$factors <- this$lApply(factor_cols, this$get_unique_values, data = data)
		names(this$factors) <- factor_cols
	}

	# Column indexes that are factors
	this$factor_indices <- which(this$col_names %in% names(this$factors))

	format_column <- function(column = stop("column: has not been specified"),
								levels = stop("levels: have not been specified")){
		as.numeric(factor(column, levels))
	}

	this$format_column <- format_column

	for(name in names(this$factors)){
		data[,name] <- this$format_column(data[,name], this$factors[[name]])
	}

	col_classes[col_classes == "character"] <- "factor"
	this$col_classes <- col_classes

	# Load UUID
	uuid <- uuid::UUIDgenerate
	this$uuid <- uuid
	this$make_name <- function(){
		paste0("_A", gsub("-", "", uuid(use.time = FALSE)))
	}

	# Alternative to make_name function does not require {uuid} package
	this$make_name_2 <- function(){
		paste0("_A", paste0(sample(c(LETTERS, 0:9), 32, TRUE), collapse = ""), collapse = "")
	}

	# TODO: Finish this
	# Default block reader
	.Reader <- function(block_name = stop("block name not specified")){
		lf <- vector(mode = "list", length = this$ncols)
		names(lf) <- this$col_names
		block <- read_df_block(this$size[block_name], block_name, lf, 
				this$ncols, this$factors, this$factor_indices - 1)
		return(block)
	}

	environment(.Reader) <- this
	this$Reader <- .Reader

	ReadBlocks1 <- function(block_names){
		lf <- vector(mode = "list", length = this$ncols)
		names(lf) <- this$col_names
		read_df_blocks(this$size[block_names], block_names, lf, 
				this$ncols, this$factors, this$factor_indices - 1)
	}
	this$ReadBlocks1 <- ReadBlocks1

	ReadBlocks <- function(block_names){
		cat("Reading ..\n")
		dfs <- this$lApply(block_names, .Reader)
		cat("Binding ..\n")
		return(do.call(rbind, dfs))
	}
	this$ReadBlocks <- ReadBlocks

	.Reader2 <- function(block_name = stop("block name not specified")){
		block <- readNumericVector(this$size[block_name], block_name)
		dim(block) <- c(length(block)/this$ncols, this$ncols)
		block <- data.frame(block)
		names(block) <- this$col_names
		for(name in names(this$factors)){
			block[,name] <- CreateFactor(block[,name], this$factors[[name]])
		}
		return(block)
	}
	this$Reader2 <- .Reader2

	.Writer <- function(block = stop("block: data to be written has not been specified"),
								block_name){
		block <- as.vector(block)
		write_numeric_vector(block, block_name)
		size <- length(block)
		names(size) <- block_name
		return(size)
	}

	environment(.Writer) <- this

	this$Writer <- .Writer

	# Function to set the file writer
	this$set_writer <- function(fun){
		environment(fun) <- this
		this$Writer <- fun
	}

	this$extn <- "bds"
	# Set file extension
	this$set_file_extension <- function(extn){
		this$extn <- extn
	}

	create_indices <- function(nrows = stop("nrows must be specified")){
		stubSize <- nrows %% this$chunkSize
	    if(stubSize != 0)
	    {
	        index <- c(1, seq(stubSize, nrows, by = this$chunkSize) + 1)
	        indices <- list()
	        for(i in 2:length(index)){
	          indices[[i-1]] <- c(index[i - 1], index[i] - 1)
	        }
	    }else{
	        index <- seq(this$chunkSize, nrows, this$chunkSize)
	        indices <- list()
	        for(i in 1:length(index))
	        {
	          indices[[i]] <- c(index[i] - this$chunkSize + 1, index[i])
	        }
	    }
	    return(indices)
	}

	environment(create_indices) <- this
	this$create_indices <- create_indices

	write_indexed_data <- function(indices = stop("indices: must be supplied"), 
		data = stop("data: must be supplied")){
		block_names <- vapply(1:length(indices), function(x)this$make_name(), "")
		block_names <- paste0(file.path(this$path, block_names), ".", this$extn)
		new_sizes <- this$lApply(1:length(block_names), function(x, data){
			this$Writer(data[indices[[x]][1]:indices[[x]][2], ], block_names[x])
		}, data = data)
		new_sizes <- unlist(new_sizes)
		this$nrows <- this$nrows + sum(new_sizes)/this$ncols
		return(new_sizes)
	}

	environment(write_indexed_data) <- this
	this$write_indexed_data <- write_indexed_data

	# Function to read a block of data
	read_block <- function(block_name){
		this$Reader(block_name)
	}
	environment(read_block) <- this
	this$read_block <- read_block
	
	this$indices <- this$create_indices(nrow(data))

	data <- as.matrix(data)
	new_sizes <- write_indexed_data(this$indices, data)
	this$size <- c(this$size, new_sizes)
	this$block_names <- c(this$block_names, names(new_sizes))

	this$append <- function(data = stop("data must be supplied")){
		if(this$isFile(data)){
			data <- this$fReader(data)
		}

		for(name in names(this$factors)){
			# Include code for appending new factor levels
			new_levels <- unique(as.character(data[,name]))
			new_levels <- new_levels[!new_levels %in% this$factors[[name]]]
			this$factors[[name]] <- c(this$factors[[name]], new_levels)
			data[,name] <- this$format_column(data[,name], this$factors[[name]])
		}

		indices <- this$create_indices(nrow(data))

		data <- as.matrix(data)
		new_sizes <- write_indexed_data(indices, data)
		this$size <- c(this$size, new_sizes)
		this$block_names <- c(this$block_names, names(new_sizes))
		save(this, file = file.path(this$path, ".meta_data"))
	}

	# Same as append but without file metadata update
	this$add <- function(data = stop("data must be supplied")){
		if(this$isFile(data)){
			data <- this$fReader(data)
		}

		for(name in names(this$factors)){
			# Include code for appending new factor levels
			new_levels <- unique(as.character(data[,name]))
			new_levels <- new_levels[!new_levels %in% this$factors[[name]]]
			this$factors[[name]] <- c(this$factors[[name]], new_levels)
			data[,name] <- this$format_column(data[,name], this$factors[[name]])
		}

		indices <- this$create_indices(nrow(data))

		data <- as.matrix(data)
		new_sizes <- write_indexed_data(indices, data)
		this$size <- c(this$size, new_sizes)
		this$block_names <- c(this$block_names, names(new_sizes))
		return(1)
	}

	this$head <- function(n = 4){
		i <- 1
		.df <- this$Reader(this$block_names[i])
		while(nrow(.df) < n && length(this$block_names) > i){
			i <- i + 1
			.df <- rbind(.df, this$Reader(this$block_names[i]))
		}
		nrowDF <- nrow(.df)
		if(n > nrowDF)
			n <- nrowDF
		.df[1:n,]
	}

	this$tail <- function(n = 4){
		i <- length(this$block_names)
		.df <- this$Reader(this$block_names[[i]])
		while(nrow(.df) < n && i > 1){
			i <- i - 1
			.df <- rbind(this$Reader(this$block_names[[i]]), .df)
		}
		.df[(nrow(.df) - n + 1):nrow(.df),]
	}

	this$n_blocks <- function(){
		length(this$block_names)
	}

	# Read the whole data set as a dataframe
	this$dataframe <- function(){
		this$ReadBlocks(this$block_names)
	}

	this$delete <- function(){
		unlink(this$path, recursive = TRUE)
		rm(list = ls(), envir = this)
		invisible()
	}

	this$save <- function(){
		save(this, file = file.path(this$path, ".meta_data"))
	}

	this$map <- function(FUN = stop("FUN: function to be applied to each block missing"), ...){
		this$lApply(this$block_names, function(block_name){
			dat <- this$Reader(block_name)
			return(FUN(dat, ...))
		})
	}

	this$copy <- function(new_path = stop("new_path: of new object must be supplied")){
		object_files <- list.files(this$path, all.files = TRUE)[-(1:2)]
		orig_files <- list.files(this$path, all.files = TRUE, full.names = TRUE)[-(1:2)]
		object_files <- file.path(new_path, object_files)
		if(!dir.exists(new_path))
			dir.create(new_path)
		file.copy(orig_files, object_files)
		new_names <- ls(this, all.names = TRUE)
		new_obj <-  lapply(new_names, get, envir = this)
		names(new_obj) <- new_names
		new_obj <- as.environment(new_obj)
		class(new_obj) <- "data_frame"
		return(new_obj)
	}

	# Creates a new data_matrix from a data frame by applying a function to it
	# Function FUN must return a matrix!
	this$new_matrix <- function(FUN = stop("FUN: transforming function not found"),
		new_path = stop("new_path: path for new matrix not found"), ...){
		matStub <- FUN(this$Reader(this$block_names[1]), ...)
		that <- data_matrix(matStub, this$chunkSize, new_path)
		new_sizes <- this$lApply(this$block_names[-1], function(block_name, ...){
			.df <- this$Reader(block_name)
			.mm <- FUN(.df, ...)
			new_name <- paste0(file.path(that$path, that$make_name()), ".", that$extn)
			# This works but should be looked at again Writer arg 1 should be a vector not a matrix
			# Unfolding occurs automatically but should it be explicit?
			return(that$Writer(.mm, new_name))
		}, ...)
		new_sizes <- unlist(new_sizes)
		that$size <- c(that$size, new_sizes)
		that$block_names <- c(that$block_names, names(new_sizes))
		this <- that
		save(this, file = file.path(this$path, ".meta_data"))
		return(this)
	}

	# Function to merge the meta data from object with the same path and data structure
	this$merge_meta <- function(object1, object2){

	}

	# Function generates new data_frame object from current object
	# This is serial for now because of factoring issues
	this$new_data_frame <- function(FUN = stop("FUN: transforming function not found"),
		new_path = stop("new_path: path for new matrix not found"), ...){
		dfStub <- FUN(this$Reader(this$block_names[1]), ...)
		that <- data_frame(dfStub, this$chunkSize, new_path)
		if(length(this$block_names > 1)){
			blockFun <- function(block_name, ...){
				.df <- FUN(this$Reader(block_name), ...)
				new_name <- paste0(file.path(that$path, that$make_name()), ".", that$extn)
				that$add(.df)
			}
			lP <- this$lApply
			this$lApply <- lapply # make this serial for now - factoring issue ...
			thats <- this$lApply(this$block_names[-1], blockFun, ...)
			this$lApply <- lP
		}
		this <- that
		save(this, file = file.path(this$path, ".meta_data"))
		return(this)
	}

	# Deletes all files not in block_names (unless hidden)
	this$cleanup <- function(){
		all_files <- list.files(this$path, full.names = TRUE)
		unlink(all_files[!(all_files %in% this$block_names)])
	}

	# Function to modify the data frame in place
	# This is serial for now because of factoring/size issues
	# Function carries out a modify map
	this$modify <- function(FUN = stop("FUN: transforming function not found"), ...){
		old_names <- this$block_names[-1]
		dfStub <- FUN(this$Reader(this$block_names[1]), ...)
		that <- data_frame(dfStub, this$chunkSize, this$path)
		if(length(old_names > 1)){
			blockFun <- function(block_name, ...){
				.df <- FUN(this$Reader(block_name), ...)
				new_name <- paste0(file.path(that$path, that$make_name()), ".", that$extn)
				that$add(.df)
			}
			lP <- this$lApply
			this$lApply <- lapply # make this serial for now - factoring issue ...
			thats <- this$lApply(old_names, blockFun, ...)
			this$lApply <- lP
		}
		unlink(old_names) # clear out old data
		this <- that
		save(this, file = file.path(this$path, ".meta_data"))
		return(this)
	}

	class(this) <- "data_frame"
	save(this, file = file.path(this$path, ".meta_data"))

	return(this)
}


#' @title print function for a data_frame
#' @description print function for a data_frame
#' @param x data_frame object to print
#' @param ... not used
#' @export
print.data_frame<- function(x, ...){
	cat(paste0("data_frame object: ncols: ", x$ncols, ", nrows: ", x$nrows, "\n"))
	print(x$head())
}

# The load_data_frame function
#' @title function to load data_frame object
#' @description function to load data_frame object
#' @param path character to folder containing object
#' @export
load_data_frame <- function(path = stop("path: to data_frame folder must be supplied")){
	this <- NULL
	load(file.path(path, ".meta_data"))
	return(this)
}
