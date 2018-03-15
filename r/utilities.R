
# Format a point for printing, based on specified precision with trailing zeros. Uniform printing for vector-like data 
# (vector, array, list).
# @param point (vector-like): nD point with floating point coordinates.
# @param precision (int): Number of digits after the decimal point.
# @return: String representation of the given point "xx.xxx yy.yyy zz.zzz...".
point2str <- function(point, precision=1)
{
    precision_str <- sprintf("%%.%df",precision)
    return(paste(lapply(point, function(x) sprintf(precision_str, x)), collapse=", "))
}


# Generate random (uniform withing bounds) nD point cloud. Dimension is based on the number of pairs in the 
# bounds input.
# @param bounds (list(vector-like)): List where each vector defines the coordinate bounds.
# @param num_points (int): Number of points to generate.
# @return (matrix): Matrix whose columns are the set of points.                         
uniform_random_points <- function(bounds, num_points)
{
    return(t(sapply(bounds, function(bnd,n=num_points) runif(n, min(bnd),max(bnd)))))
}
                                 

# Distances between points transformed by the given transformation and their
# location in another coordinate system. When the points are only used to evaluate
# registration accuracy (not used in the registration) this is the target registration
# error (TRE).
# @param tx (SimpleITK transformation): Transformation applied to the points in point_data
# @param point_data (matrix): Matrix whose columns are points which we transform using tx.
# @param reference_point_data (matrix): Matrix whose columns are points to which we compare 
#                                       the transformed point data.                                              
# @return (vector): Distances between the transformed points and the reference points.
target_registration_errors <- function(tx, point_data, reference_point_data)
{
    transformed_points_mat <- apply(point_data, MARGIN=2, tx$TransformPoint)
    return (sqrt(colSums((transformed_points_mat - reference_point_data)^2)))
}
                                 
                                 
# Check whether two transformations are "equivalent" in an arbitrary spatial region 
# either 3D or 2D, [x=(-10,10), y=(-100,100), z=(-1000,1000)]. This is just a sanity check, 
# as we are just looking at the effect of the transformations on a random set of points in
# the region.
print_transformation_differences <- function(tx1, tx2)
{
    if (tx1$GetDimension()==2 && tx2$GetDimension()==2)
    {
        bounds <- list(c(-10,10), c(-100,100))
    }
    else if(tx1$GetDimension()==3 && tx2$GetDimension()==3)
    {
        bounds <- list(c(-10,10), c(-100,100), c(-1000,1000))
    }
    else
        stop('Transformation dimensions mismatch, or unsupported transformation dimensionality')
    num_points <- 10
    point_data <- uniform_random_points(bounds, num_points)
    tx1_point_data <- apply(point_data, MARGIN=2, tx1$TransformPoint)
    differences <- target_registration_errors(tx2, point_data, tx1_point_data)
    cat(tx1$GetName(), "-", tx2$GetName(), ":\tminDifference: ", 
        toString(min(differences)), " maxDifference: ",toString(max(differences))) 
}

#
# This function displays the effects of the deformable transformation on a grid of points by scaling the
# initial displacements (either of control points for BSpline or the deformation field itself). It does
# assume that all points are contained in the range(-2.5,-2.5), (2.5,2.5) - for display.
#
display_displacement_scaling_effect <- function(s, original_x_mat, original_y_mat, tx, original_control_point_displacements)
{
    if(tx$GetDimension()!=2)
        stop('display_displacement_scaling_effect only works in 2D')

    tx$SetParameters(s*original_control_point_displacements)
    transformed_points <- mapply(function(x,y) tx$TransformPoint(c(x,y)), original_x_mat, original_y_mat)
        
    plot(original_x_mat,original_y_mat, xlim=c(-2.5,2.5), ylim=c(-2.5,2.5), pch=19, col="blue", xlab="", ylab="", las=1)
    points(transformed_points[1,], transformed_points[2,], col="red", pch=17)
    legend('top', col= c("red", "blue"), pch=c(17,19), legend = c("transformed points", "original points"))
}

   
# Create a list representing a regular sampling of the parameter space.
# Args:two or more vectors representing parameter values. The order
#      of the vectors should match the ordering of the SimpleITK transformation
#      parameterization (e.g. Similarity2DTransform: scaling, rotation, tx, ty)
# Return:
# List of  vectors representing the regular grid sampling.
parameter_space_regular_grid_sampling <- function(...) {
    df <- expand.grid(list(...))
    return(lapply(split(df,seq_along(df[,1])), as.vector))
}


# Create a list representing a regular sampling of the 3D similarity transformation parameter space. As the
# SimpleITK rotation parameterization uses the vector portion of a versor we don't have an
# intuitive way of specifying rotations. We therefor use the ZYX Euler angle parametrization and convert to
# versor.
# Args:
#   thetaX, thetaY, thetaZ: vectors with the Euler angle values to use.
#   tx, ty, tz: vectors with the translation values to use.
#   scale: vector with the scale values to use.
#    Return:
#        List of vectors representing the parameter space sampling (vx,vy,vz,tx,ty,tz,s).
similarity3D_parameter_space_regular_sampling <- function(thetaX, thetaY, thetaZ, tx, ty, tz, scale) {
  euler_sampling <- parameter_space_regular_grid_sampling(thetaX, thetaY, thetaZ, tx, ty, tz, scale)
  # replace Euler angles with quaternion vector
  for(i in seq_along(euler_sampling)) {
    euler_sampling[[i]][1:3] <- eul2quat(as.double(euler_sampling[[i]][1]),
                                         as.double(euler_sampling[[i]][2]),
                                         as.double(euler_sampling[[i]][3]))
  }
  return(euler_sampling)
}

# Translate between Euler angle (ZYX) order and quaternion representation of a rotation.
# Args:
#   ax: X rotation angle in radians.
#   ay: Y rotation angle in radians.
#   az: Z rotation angle in radians.
#   atol: tolerance used for stable quaternion computation (qs==0 within this tolerance).
# Return:
#      Vector with three entries representing the vectorial component of the quaternion.
eul2quat <- function(ax, ay, az, atol=1e-8) {
    # Create rotation matrix using ZYX Euler angles and then compute quaternion using entries.
    cx <- cos(ax)
    cy <- cos(ay)
    cz <- cos(az)
    sx <- sin(ax)
    sy <- sin(ay)
    sz <- sin(az)
    r <- array(0,c(3,3))
    r[1,1] <- cz*cy
    r[1,2] <- cz*sy*sx - sz*cx
    r[1,3] <- cz*sy*cx+sz*sx

    r[2,1] <- sz*cy
    r[2,2] <- sz*sy*sx + cz*cx
    r[2,3] <- sz*sy*cx - cz*sx

    r[3,1] <- -sy
    r[3,2] <- cy*sx
    r[3,3] <- cy*cx

    # Compute quaternion:
    qs <- 0.5*sqrt(r[1,1] + r[2,2] + r[3,3] + 1)
    qv <- c(0,0,0)
    # If the scalar component of the quaternion is close to zero, we
    # compute the vector part using a numerically stable approach
    if(abs(qs - 0.0) < atol) {
        i <- which.max(c(r[1,1], r[2,2], r[3,3]))
        j <- (i+1)%%3 + 1
        k <- (j+1)%%3 + 1
        w <- sqrt(r[i,i] - r[j,j] - r[k,k] + 1)
        qv[i] <- 0.5*w
        qv[j] <- (r[i,j] + r[j,i])/(2*w)
        qv[k] <- (r[i,k] + r[k,i])/(2*w)
    } else {
        denom <- 4*qs
        qv[1] <- (r[3,2] - r[2,3])/denom;
        qv[2] <- (r[1,3] - r[3,1])/denom;
        qv[3] <- (r[2,1] - r[1,2])/denom;
    }
    return(qv)
}
