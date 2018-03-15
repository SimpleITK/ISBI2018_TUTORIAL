library(ggplot2)

# Labels used in the POPI dataset
popi_body_label <- 0
popi_air_label <- 1
popi_lung_label <- 2


# Callback invoked when the StartEvent happens, sets up our new data.
# Functions ending in _jn are for use with Jupyter notebooks, as the display
# behaviour is a bit different.
# Note that we could use a special plotting environment instead of the global
# environment, but we may want to use the metrics elsewhere and they are easier to
# access if they are global
start_plot <- function()
{   #global empty vectors (via assignment operator) 
    metric_values <<- c()
    multires_iterations <<- c()
}


end_plot_jn <- function()
{    
    Multi <- rep(NA, length(metric_values))
    Multi[multires_iterations] <- "M"
    DDF <- data.frame(IterationNumber=1:length(metric_values), 
                      MetricValue=metric_values,
                      MultiresIteration=Multi)
    DDFM <- subset(DDF, !is.na(MultiresIteration))
    pl <- ggplot(DDF, aes(x=IterationNumber, y=MetricValue)) + 
          geom_line() + 
          geom_point(data=DDFM, aes(colour=MultiresIteration)) +
          theme(legend.position="none")    
    print(pl)
    rm(metric_values, pos = ".GlobalEnv")
    rm(multires_iterations, pos = ".GlobalEnv")
}


# Callback invoked when the IterationEvent happens, update our data and display new figure.
# Note that this won't appear as an animation in R studio, but you can use the arrows to cycle
# through plots
plot_values <- function(registration_method)
{
    metric_values <<- c(metric_values, registration_method$GetMetricValue())
    Multi <- rep(NA, length(metric_values))
    Multi[multires_iterations] <- "M"
    DDF <- data.frame(IterationNumber=1:length(metric_values),
                      MetricValue=metric_values,
                      MultiresIteration=Multi)
    DDFM <- subset(DDF, !is.na(MultiresIteration))

    pl <- ggplot(DDF, aes(x=IterationNumber, y=MetricValue)) +
          geom_line() +
          theme(legend.position="none")

    if(nrow(DDFM) > 1) {
        pl <- pl + geom_point(data=DDFM, aes(colour=MultiresIteration))
    }
    print(pl)
    dev.flush()
    Sys.sleep(0)
}


# Use this one inside a notebook
plot_values_jn <- function(registration_method)
{    
    ## No point attempting to plot every one in a notebook
    metric_values <<- c(metric_values, registration_method$GetMetricValue())
}


# Callback invoked when the sitkMultiResolutionIterationEvent happens, update the index into the 
# metric_values list. 
update_multires_iterations <- function()
{
    multires_iterations <<- c(multires_iterations, length(metric_values)+1)
}


#
# Get a coronal slice with overlaid contour of the mask for the specific slice index in all temporal images.
#
temporal_coronal_with_overlay <- function(coronal_slice_index, images, masks, label, window_min, window_max)
{
    # Extract the 2D images and masks.
    slices <- lapply(images, function(img, slc) img[,slc,], slc=coronal_slice_index)
    slice_masks <- lapply(masks, function(msk, slc, lbl) msk[,slc,]==lbl , slc=coronal_slice_index, lbl=label)

    # Resample the image (linear interpolation) and mask (nearest neighbor interpolation) into an isotropic grid,
    # required for display.
    original_spacing <- slices[[1]]$GetSpacing()
    original_size <- slices[[1]]$GetSize()
    min_spacing <- min(original_spacing)
    new_spacing <- c(min_spacing, min_spacing)
    new_size <- c(as.integer(round(original_size[1]*(original_spacing[1]/min_spacing))),
                  as.integer(round(original_size[2]*(original_spacing[2]/min_spacing))))
    resampled_slices <- lapply(slices, function(slc, sz, spc) Resample(slc, sz, Transform(),
                                                                       "sitkLinear", slc$GetOrigin(),
                                                                       spc, slc$GetDirection(), 0.0,
                                                                       slc$GetPixelID()), sz=new_size, spc=new_spacing)
    resampled_slice_masks <- lapply(slice_masks, function(msk, sz, spc) Resample(msk, sz, Transform(),
                                                                       "sitkNearestNeighbor", msk$GetOrigin(),
                                                                       spc, msk$GetDirection(), 0.0,
                                                                       msk$GetPixelID()), sz=new_size, spc=new_spacing)

    # Create the overlay: cast the mask to expected label pixel type, and do the same for the image after
    # window-level, accounting for the high dynamic range of the CT.
    overlaid_slices <- mapply( function(slc, msk, win_min, win_max) LabelMapContourOverlay(Cast(msk, "sitkLabelUInt8"),
                                                                    Cast(IntensityWindowing(slc,
                                                                                            windowMinimum=win_min,
                                                                                            windowMaximum=win_max),
                                                                         "sitkUInt8"),
                                                                     opacity = 1,
                                                                     c(0,0), c(2,2)),
                                resampled_slices,
                                resampled_slice_masks, win_min=window_min, win_max=window_max)
    # Create the temporal slice, 3D volume representing 2D coronal+time
    temporal_image <- Image(c(overlaid_slices[[1]]$GetSize(), length(overlaid_slices)), overlaid_slices[[1]]$GetPixelID())
    # Two subtle points: (1) to paste the 2D slice into the 3D volume we need to make it a 3D slice (JoinSeries),
    #                    (2) the Paste function uses SimpleITK indexing, requiring the seq()-1.
    invisible(mapply(function(slice, index) temporal_image<<- Paste(temporal_image, JoinSeries(slice), c(slice$GetSize(),1), c(0,0,0), c(0,0,index)),
                     overlaid_slices, seq(length(overlaid_slices))-1))

    return(temporal_image)
}
