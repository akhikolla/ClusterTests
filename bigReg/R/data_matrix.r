#' Creates a data_matrix object
#' @title function to create a data_frame object
#' @description function to create a data_matrix object. The data_matrix object
#'              is an object that is held on disk. It is written to a folder
#'              \code{path} on disk where the data is written to in blocks or
#'              chunks. The data is written in binary format using a C++ function
#'              in purely numerical data.
#' @param data object to be converted into a data_matrix object
#' @param chunkSize number of rows to be used in each chunk
#' @param path character to folder where the object will be created
#' @param nCores the number of cores to use defaults to parallel::detectCores()
#' @param ... not used at the moment
#' @export
data_matrix <- function(data = stop("data: matrix must be supplied"),
	chunkSize = stop("chunkSize must be specified, a good number is 50000"),
	path = stop("path must be specified"), nCores = parallel::detectCores(), ...){

	this <- new.env()
	this$size <- vector(mode = "numeric", length = 0)

	if(!dir.exists(path)){
		dir.create(path)
	}

	this$path <- path

	# Default lapply function
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

	set_reader <- function(fun){
		environment(fun) <- this
		this$Reader <- fun
	}
	environment(set_reader) <- this
	this$set_reader <- set_reader

	# Common file reader
	this$fReader <- read.csv

	set_file_reader <- function(fun){
		environment(fun) <- this
		this$fReader <- fun
	}

	environment(set_file_reader) <- this
	this$set_file_reader <- set_file_reader

	# This reads the file if it is a data path
	if(this$isFile(data)){
		data <- as.matrix(this$fReader(data))
	}

	this$ncols <- ncol(data)
	this$col_names <- dimnames(data)[[2]]

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
		block <- read_matrix_block(block_name, this$size[block_name], this$ncols)
		dimnames(block)[[2]] <- this$col_names
		return(block)
	}

	environment(.Reader) <- this
	this$Reader <- .Reader

	ReadBlocks1 <- function(block_names){
		blocks <- read_matrix_blocks(block_names, this$size[block_names], this$ncols)
		dimnames(blocks)[[2]] <- this$col_names
		return(blocks)
	}

	this$ReadBlocks1 <- ReadBlocks1

	ReadBlocks <- function(block_names){
		cat("Reading ..\n")
		mats <- this$lApply(block_names, .Reader)
		cat("Binding ..\n")
		return(do.call(rbind, mats))
	}
	this$ReadBlocks <- ReadBlocks

	ReadAsList <- function(){
		mats <- this$lApply(this$block_names, .Reader)
	}
	this$ReadAsList <- ReadAsList

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

	# TODO: Something else other than a for loop?
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
		return(new_sizes)
	}

	environment(write_indexed_data) <- this
	this$write_indexed_data <- write_indexed_data
	
	this$indices <- this$create_indices(nrow(data))

	new_sizes <- write_indexed_data(this$indices, data)
	this$size <- c(this$size, new_sizes)
	this$block_names <- c(this$block_names, names(new_sizes))

	this$append <- function(data = stop("data: matrix must be supplied")){
		if(this$isFile(data)){
			data <- as.matrix(this$fReader(data))
		}

		indices <- this$create_indices(nrow(data))

		data <- as.matrix(data)
		new_sizes <- write_indexed_data(indices, data)
		this$size <- c(this$size, new_sizes)
		this$block_names <- c(this$block_names, names(new_sizes))
		save(this, file = file.path(this$path, ".meta_data"))
	}

	# Same as append but without file metadata update
	this$add <- function(data = stop("data: matrix must be supplied")){
		if(this$isFile(data)){
			data <- as.matrix(this$fReader(data))
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
		nrowDF <- nrow(.df)
		if(n > nrowDF)
			n <- nrowDF
		.df[(nrow(.df) - n + 1):nrow(.df),]
	}

	this$n_blocks <- function(){
		length(this$block_names)
	}

	# Read the whole data set as a matrix
	this$matrix <- function(){
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
	this$new_matrix <- function(FUN = stop("FUN: transforming function not found"),
		new_path = stop("new_path: path for new matrix not found"), ...){
		matStub <- FUN(this$Reader(this$block_names[1]), ...)
		that <- data_matrix(matStub, this$chunkSize, new_path)
		if(length(this$block_names) > 1){
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
		}
		save(that, file = file.path(that$path, ".meta_data"))
		return(that)
	}

	# Function generates new data_frame object from current object
	this$new_data_frame <- function(FUN = stop("FUN: transforming function not found"),
		new_path = stop("new_path: path for new matrix not found"), ...){
		dfStub <- FUN(this$Reader(this$block_names[1]), ...)
		that <- data_frame(dfStub, this$chunkSize, new_path)
		if(length(this$block_names) > 1){
			blockFun <- function(block_name, ...){
				.df <- FUN(this$Reader(block_name), ...)
				new_name <- paste0(file.path(that$path, that$make_name()), ".", that$extn)
				that$add(.df)
			}
			lP <- this$lApply
			this$lApply <- lapply # make this serial for now - factoring issue ...
			this$lApply(this$block_names[-1], blockFun, ...)
			this$lApply <- lP
		}
		save(that, file = file.path(that$path, ".meta_data"))
		return(that)
	}

	# Deletes all files not in block_names (unless hidden)
	this$cleanup <- function(){
		all_files <- list.files(this$path, full.names = TRUE)
		unlink(all_files[!(all_files %in% this$block_names)])
	}

	# Function to modify the data frame in place
	# This is serial for now because of factoring/size issues
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
		save(that, file = file.path(that$path, ".meta_data"))
		return(that)
	}

	class(this) <- "data_matrix"

	save(this, file = file.path(this$path, ".meta_data"))

	return(this)
}


#' @title print function for a data_matrix
#' @description print function for a data_matrix
#' @param x data_matrix object to print
#' @param ... not used
#' @export
print.data_matrix<- function(x, ...){
	cat(paste0("data_matrix object: ncols: ", x$ncols, ", nrows: ", x$nrows, "\n"))
	print(x$head())
}


# The load_data_frame function
#' @title function to load data_frame object
#' @description function to load data_frame object
#' @param path character to folder containing object
#' @export
load_data_matrix <- function(path = stop("path: to data_matrix folder must be supplied")){
	this <- NULL
	load(file.path(path, ".meta_data"))
	return(this)
}
