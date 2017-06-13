#!/usr/bin/env octave

## Copyright (C) 2014-2017 David Pinto <david.pinto@bioch.ox.ac.uk>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.

pkg load image;
pkg load optim;
pkg load statistics;

pkg load bioformats;
javaMethod ("enableLogging", "loci.common.DebugTools", "ERROR");

page_output_immediately (1);
more ("off");


function [img, voxel_sizes] = read_image (fpath)
  reader = bfGetReader (fpath);
  n_planes = reader.getImageCount ();
  img = bfGetPlane (reader, 1);
  img = postpad (img, n_planes, 0, 4);
  for p_idx = 2:n_planes
    img(:,:,p_idx) = bfGetPlane (reader, p_idx);
  endfor

  ## Dimensions in the 3rd dimension are in the order [channel, z, time]
  ## but we don't have, nor are interested, in channel.
  img = reshape (img, [size(img)(1:2) reader.getSizeZ() reader.getSizeT()]);

  metadata_store = reader.getMetadataStore ();
  micron = java_get ("ome.units.UNITS", "MICROM");
  voxel_sizes = [
    metadata_store.getPixelsPhysicalSizeX(0).value(micron)
    metadata_store.getPixelsPhysicalSizeY(0).value(micron)
    metadata_store.getPixelsPhysicalSizeZ(0).value(micron)
  ]';
endfunction

## Compute threshold value using triangle method adapted for
##fluorescence microscopy.
##
## The algorithm was designed for "normal" light microscopy images
##where the histogram has a high peak for pixels of white intensity
##and the part of interest is the long tail after the peak.
##
## Our histogram is the other way around (note that on the publication
## they put white on the left though), which is the adaptation that we
## do.
##
## See discussion of the algorithm at
##http://forum.imagej.net/t/understanding-imagej-implementation-of-the-triangle-algorithm-for-threshold/752/7
##
## TODO this should be done in the Image package
##
## FIXME this does not return the same threshold values as ImageJ.
##       Not sure what is causing the difference though.  But the
##       different threshold values seem to work even better for our
##       images.
function [thresh] = triangle_threshold (img)

  nbins = 256;
  [counts, bins_centers] = hist (img(:), nbins);

  ## The algorithm assumes a high peak for bright values but this is
  ## fluorescence microscopy so invert the histogram.
  counts = counts(end:-1:1);
  [h_max, h_max_idx] = max (counts);

  ## The reference for the algorithm does not specify if the min is
  ## the first non-0 value or the first 0 zero before it.  But ImageJ
  ## and "Fundamentals of Image Processing" say bmin = (p=0)% so we
  ## "draw" the triangle from the non-existing bar before the first
  ## bar, to the highest histogram peak.  And we use the center of
  ## each bar for computations.
  triangle_range = 1:(h_max_idx-1);
  norm_heights = counts(triangle_range) ./ h_max;

  ## The positions of each bin center, as distance from the center of
  ## the 0 height bin before the histogram.  Do not compute the last
  ## bin because the threshold won't be there.
  norm_bin_positions = (1 ./ h_max_idx) .* triangle_range;

  distances = norm_bin_positions - norm_heights;
  [~, thresh] = max (distances);

  ## The space to find the threshold has been reduced twice. First, we
  ## throw away the whole dynamic range and use only the part where
  ## there's counts.  Then from that histogram we use only the part of
  ## the histogram from the peak to one of the sides.  Now we undo all
  ## that to get a threshold in range [0 1].
  range_end = getrangefromclass (img)(2);
  thresh_offset = bins_centers(1) ./ range_end;
  hist_range = (bins_centers(end) - bins_centers(1)) ./ range_end;

  hist_thresh = (numel (counts) - thresh) ./ hist_range;
  thresh = thresh_offset + hist_thresh;
  thresh = im2double (cast (thresh, class (img)));
endfunction

## Finds the initial coordinates for centriole centers by finding the
## weighted centroids of the image peaks.
function [obj_cen] = weighted_centroids (img)
  if (! isinteger (img))
    ## we could just use abs() but if data is floating point, we may
    ## be sure we are not messing up something.
    error ("Data is not integer class");
  endif

  conn = true (3, 3, 3);

  ## Mask for each centriole in image.
  mask = im2bw (img, triangle_threshold (img));
  mask = bwareaopen (mask, 20, conn);

  ## Split touching centrioles using watershed lines.
  img(! mask) = 0;
  mask(! watershed (imcomplement (img), conn)) = false;

  props = regionprops (bwconncomp (mask, conn), img, {"WeightedCentroid"});
  obj_cen = cell2mat ({props(:).WeightedCentroid}(:));
  obj_cen(:,[1 2]) = obj_cen(:,[2 1]); # swap x,y to get row,column
endfunction

## Create a structuring element ball shaped, for elements below the
## cutoff distance.
##
## I wanted to do 'fspecial ("disk", [x y z]) > 0' but fspecial only
## works up to 2 dimensions and is not trivial at the moment to extend
## it on the Image package because of the way the border pixels need
## to be handled there.
function [se] = create_ball_se (voxel_sizes, cutoff)
  cutoff_voxels = round (cutoff ./ voxel_sizes);
  limits = voxel_sizes .* cutoff_voxels;
  nd = numel (voxel_sizes);
  X = cell (nd, 1);
  for i = 1:nd
    X{i} = linspace (-limits(i), limits(i), cutoff_voxels(i)*2 +1);
  endfor
  ## The cell array in this shape allows for the cell2mat afterwards.
  grids = cell ([ones(1, nd) nd]);
  [grids{:}] = ndgrid (X{:});
  grids = cell2mat (grids);
  distances = sqrt (sum (grids.^2, 4));
  se = distances < cutoff;
endfunction

## Group centroids in distances less than cutoff_distance.  Returns a
## cell array, each element a group of rows from points within the
## cutoff distance.
##
## We should be using pdist, linkage, and cluster from the statistics
## package:
##
##    d = pdist (p);
##    z = linkage (d);
##    cidx = cluster (z, "cutoff", 0.5, "criterion", "distance");
##
## but the cluster function is still not implemented.  So we fudge it
## with mathematical morphology.
function [groups, n_labels] = cluster_points (points, sz, voxel_sizes,
                                              cutoff_distance)
  canvas = false (sz);
  points_linear = sub2ind (sz, num2cell (round (points), 1){:});
  canvas(points_linear) = true;
  se = create_ball_se (voxel_sizes, cutoff_distance);
  canvas = imdilate (canvas, se);
  [labelled, n_labels] = bwlabeln (canvas);
  groups = labelled(points_linear);
  groups = accumarray (groups, num2cell (points, 2),
                       [n_labels 1], @(x) {x});
endfunction

## Generate coordinates for pixels in physical space.
##
## LENGTHS is a vector of integers specifying the length of each
##   dimension in pixels.
##
## RESOLUTION is a vector specifying the length of each element on
##   each dimension.
##
## Returns a matrix of size KxD where K is the number of voxels and D
## is the number of dimensions.
function [coords] = get_coordinates (lengths, resolution)
  nd = numel (lengths);
  n_voxels = prod (lengths);

  strides = cell (nd, 1);
  for di = 1:nd
    strides{di} = linspace (0, resolution(di).*lengths(di), lengths(di));
  endfor

  coords = cell (1, nd);
  [coords{:}] = ndgrid (strides{:});
  coords = cell2mat (cellfun (@vec, coords, 'UniformOutput', false));
endfunction

## Create image with gaussian functions.
##
## LENGTHS is a vector of integers specifying the output size.
##
## SIGMAS is the sigma for each gaussian.  We assume the same value
##   for all dimensions (a blob).
##
## CENTERS is KxN arrays where K is the number of gaussians and N is
##   the number of dimensions (numel of LENGTHS).  Coordinates are
##   from the top left corner and physical coordinates (not indices).
##
## HEIGHTS is a vector with the peak height of each gaussian.
##
## RESOLUTION is a vector specifying the length of each element on
##   each dimension.
function [gaussian] = draw_gaussians (lengths, centers, sigmas,
                                      heights, resolution)
  n_gaussians = rows (centers);
  nd = numel (lengths);

  init_dist = cell (nd, 1);
  for di = 1:nd
    init_dist{di} = linspace (0, resolution(di).*lengths(di), lengths(di));
  endfor

  dists = cell (n_gaussians, 1);
  for gi = 1:n_gaussians
    dist = 0;
    for di = 1:nd
      dist = dist .+ (vec (init_dist{di} - centers(gi,di), di) .^2);
    endfor
    dists{gi} = dist;
  endfor

  gaussian = zeros (lengths);
  for gi = 1:n_gaussians
    gaussian += heights(gi) .* exp (- (dists{gi} ./ (2.*(sigmas(gi).^2))));
  endfor
endfunction

## Add multiple gaussians.
##
## PARAMS is an array with the gaussian values.  It should be a
##   (2+ND)xN matrix, one column for each gaussian.  Each column has
##   the paramaters for a gaussian in the order:
##
##       [HEIGHT; CENTERS; SIGMA]
##
##   It assumes the same sigma for all dimensions while CENTERS is a
##   NDx1 array.
##
## X is a matrix of size KxND, where each row corresponds to the
##   coordinates where to evaluate the function.  K is then the number
##   of points to evaluate the function, and ND is the number of
##   dimensions.
##
## N is the number of gaussians to use.  It also defines the expected
##   numel of PARAMS.
##
## Background is assumed to be zero so add whatever to Y.

function [y] = gaussians (params, x, n = 1)
  nd = columns (x);

  ## Each column are the parameters for one more gaussian.
  if (any (size (params) != [2+nd n]))
    error ("gaussians: wrong sizes");
  endif

  heights = params(1,:);
  centers = params(2:end-1,:);
  sigmas = params(end,:);

  ## A single 3D gaussian is of the form:
  ##
  ##          (  (  (x-xo)^2 )   (  (y-yo)^2 )   (  (z-zo)^2 ))
  ##    A*exp (- (-----------) + (-----------) + (-----------))
  ##          (  (  2*sx^2   )   (  2*sy^2   )   (  2*sz^2   ))
  ##
  ## But because we use the same sigma for all dimensions, we can
  ## instead have:
  ##
  ##          (  (  (x-xo)^2 + (y-yo)^2 + (z-zo)^2  ))
  ##    A*exp (- (----------------------------------))
  ##          (  (              2*s^2               ))
  ##

  y = zeros (rows (x), 1);
  for i = 1:n
    dists = x - (centers(:,i).');
    y += heights(i) .* exp (- (sum (dists.^2, 2) ./ (2 .* sigmas(i).^ 2)));
  endfor
endfunction

function [centres, fitted] = fit_to_gaussians (data, guess, voxel_sizes)

  data = double (data);

  ## Initial guesses.
  peaks = data(sub2ind (size (data), num2cell (round (guess), 1){:}));
  sigma_guess = 0.05;

  ## To physical dimensions
  guess -= 1;
  guess .*= voxel_sizes;

  x0 = [peaks(1)       peaks(2)
        guess(1,:)(:)  guess(2,:)(:)
        sigma_guess    sigma_guess];

  lower_bounds = [x0(1,:)*0.9 # peak is off by 10%
                  x0([2 3 4],:)-0.1 # coordinates off by 0.1µm
                  x0(5,:)-0.1]; # sigma off by 0.1

  upper_bounds = [x0(1,:)*1.1 # peak is off by 10%
                  x0([2 3 4],:)+0.1 # coordinates off by 0.1µm
                  x0(5,:)+0.1]; # sigma off by 0.1

  avg_kernel = fspecial ("average", [5 5]);
  background = min (convn (data, avg_kernel, "valid")(:));

  ## Function to fit two 3d gaussians.
  g2d = @(p, x) gaussians (reshape (p, [5 2]), x, 2) + background;

  ## Not really x, one column per dimension, one row per voxel.
  xdata = get_coordinates (size (data), voxel_sizes);
  ydata = data(:);

  [x, ~, r, flag] = lsqcurvefit (g2d, x0(:), xdata, ydata, lower_bounds(:),
                                 upper_bounds(:));

  fitted = reshape (uint16 (g2d (x, xdata)), size (data));
  centres = reshape (x, [5 2])(2:4,:).';

  ## Back to pixel coordinates.
  centres ./= voxel_sizes;
  centres += 1;
endfunction

##
##
##

# Note that centroids coordinates are [x, y, z]
function m = centroids_mask (centroids, dims)
  m   = false (dims);
  cr  = round (centroids);
  ## Centroids are [x y ...] and we need [row cols ...]
  ind = sub2ind (dims, cr(:,2), cr(:,1), num2cell (cr(:,3:end), 1){:});
  m(ind) = true;
endfunction

## Create an array for writing with varargin interleaved.
##
## Each array on varargin should be of size MxNxK.
function imwrite_log (fpath, varargin)
  mlog = cat (3, cellfun (@im2uint16, varargin, "UniformOutput", false){:});
  mlog_size = size (mlog);
  mlog = reshape (mlog, [mlog_size(1:2) 1 mlog_size(3)]);
  imwrite (mlog, fpath);
endfunction

function [status] = main (fpath, log_fpath)
  narginchk (1, 2);
  cutoff_distance = 0.5; # distance in microns

  [img, voxel_sizes] = read_image (fpath);
  cutoff_voxels = round (cutoff_distance ./ voxel_sizes);
  sz = size (img);
  nd = ndims (img);

  initial_coords = weighted_centroids (img);
  [coords_groups, n_groups] = cluster_points (initial_coords, sz, voxel_sizes,
                                              cutoff_distance);

  for i = 1:n_groups
    centroids = cell2mat (coords_groups{i});
    n_centroids = rows (centroids);
    if (n_centroids < 2)
      continue; # nothing yet
    elseif (n_centroids == 2)
      min_d = min (centroids, [], 1);
      max_d = max (centroids, [], 1);
      init = max (round (min_d - cutoff_voxels), 1);
      ends = min (round (max_d + cutoff_voxels), sz);
      idx = cell (nd, 1);
      for dim = 1:nd
        idx{dim} = init(dim):ends(dim);
      endfor
      snippet = img(idx{:});
    else # > 2
      continue; # nothing yet
    endif

    guess = centroids - init + 1;
    [centres, fit] = fit_to_gaussians (snippet, guess, voxel_sizes);

    fname = sprintf ("%s%i.tif", fpath, i);
    imwrite_log (fname, snippet, fit);
    ## mark = max (snippet(:));
    ## mark += double (mark) * 1.25;
    ## snippet(sub2ind (size (snippet), num2cell (round (centres), 1){:})) = mark;
    ## snippet = reshape (snippet, rows (snippet), columns (snippet), 1,
    ##                    size (snippet, 3));
    ## imwrite (snippet, name);
  endfor


  ## [mask, centroids] = get_centroids (img);

  ## if (nargin > 1) # only create log image if requested
  ##   imwrite_log (log_fpath, imadjust (img), mask,
  ##                centroids_mask (centroids, size (img)));
  ## endif

  ## centroids(:,1:3) .*= voxel_sizes; # to physical dimensions

  ## printf ("x, y, z, t\n");
  ## printf ("%f,%f,%f, %i\n", centroids');

  status = 0;
endfunction

main (argv (){:});
