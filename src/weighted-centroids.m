#!/usr/bin/env octave

## Copyright (C) 2014-2016 David Miguel Susano Pinto
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

## TODO: When we implement yen on Octave, we can do the mask ourselves
##       instead of relying on ImageJ.
pkg load imagej;

pkg load bioformats;
javaMethod ('enableLogging', 'loci.common.DebugTools', 'ERROR');


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


function [mask, obj_cen] = get_centroids (img)
  if (! isinteger (img))
    ## we could just use abs() but if data is floating point, we may
    ## be sure we are not messing up something.
    error ("data is not of integer type");
  endif

  ## Connectivity.  Analyze time-series but do not connect objects over time.
  conn = true (3, 3, 3);

  ## Mask for each centriole in image.
  mask = im2bw (img, ij_threshold (img, "Yen"));
  mask = bwareaopen (mask, 20, conn);

  ## Split touching centrioles using watershed lines.
  img(! mask) = 0;
  mask(! watershed (imcomplement (img), conn)) = false;

  props = regionprops (bwconncomp (mask, conn), img, {"WeightedCentroid"});
  obj_cen = cell2mat ({props(:).WeightedCentroid}(:));
endfunction

## Note that centroids coordinates are [x, y, z]
function m = centroids_mask (centroids, dims)
  m   = false (dims);
  cr  = round (centroids);
  ## Centroids are [x y ...] and we need [row cols ...]
  ind = sub2ind (dims, cr(:,2), cr(:,1), num2cell (cr(:,3:end), 1){:});
  m(ind) = true;
endfunction

function imwrite_log (fpath, varargin)
  mlog = cat (3, cellfun (@im2uint8, varargin, "UniformOutput", false){:});
  mlog_size = size (mlog);
  mlog = reshape (mlog, [mlog_size(1:2) 1 prod(mlog_size(3:end))]);
  imwrite (mlog, fpath);
endfunction


function [status] = main (fpath_img, fpath_log)
  narginchk (1, 2);

  [img, voxel_sizes] = read_image (fpath_img);
  [mask, centroids] = get_centroids (img);

  if (nargin > 1) # only create log image if requested
    imwrite_log (fpath_log, imadjust (img, stretchlim (img(:), 0)), mask,
                 centroids_mask (centroids, size (img)));
  endif

  centroids(:,1:3) .*= voxel_sizes; # to physical dimensions

  printf ("x\ty\tz\n");
  printf ("%f\t%f\t%f\n", centroids');

  status = 0;
endfunction


main (argv (){:});
