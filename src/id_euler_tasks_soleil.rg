import "regent"

local c = regentlib.c
local std = terralib.includec("stdlib.h")
local sqrt = regentlib.sqrt(double)
local abs = regentlib.fabs(double)
local ceil = regentlib.ceil(double)
local fmod = regentlib.fmod(double)

--local blas = require("BLAS_functions")
local DBL_EPSILON = 2.220446e-16

return function(qoi_fld, field_space)

local id_exp = {}

--task to find absolute tolerance for MGSQR
task id_exp.get_abs_tol(domain : region(ispace(int3d), field_space), tolerance : region(ispace(int2d), double), rel_tol : double)
where reads(domain.[qoi_fld]), writes(tolerance)
do

  var bounds = tolerance.bounds.lo

  --find maximum value in the domain
  var sum_val : double = 0.0

  for i in domain do
    sum_val += abs(domain[i].[qoi_fld])
  end 

  tolerance[bounds] = sum_val * rel_tol

end

--task to flatten a 3d region to a 1d vector for subsampled algorithms
task id_exp.flatten_3d_to_vector_subsampled(region_3d : region(ispace(int3d), field_space), vec_region : region(ispace(int2d), double), 
                                            column_id : int, sub_ind_x1 : int, sub_ind_x2 : int, sub_ind_x3 : int)
where reads(region_3d.[qoi_fld]), writes(vec_region) 
do

  var bounds = region_3d.bounds
  var entry_x : int = vec_region.bounds.lo.x
  var entry_y : int  = column_id

  for i = bounds.lo.x, bounds.hi.x + 1, sub_ind_x1 do
    for j = bounds.lo.y, bounds.hi.y + 1, sub_ind_x2 do
      for k = bounds.lo.z, bounds.hi.z + 1, sub_ind_x3 do
        vec_region[{x = entry_x, y = entry_y}] = region_3d[{x = i, y = j, z = k}].[qoi_fld]
        entry_x += 1
      end
    end
  end

end

--task to perform trilinear interpolation on a subsampled 3d region
task id_exp.trilinear_interp(A_c : region(ispace(int2d), double), S : region(ispace(int2d), double), k : region(ispace(int2d), int), 
                             nx : int, ny : int, nz : int, dx : double, dy : double, dz : double, nx_total : int, ny_total : int, nz_total : int, subsampling_x : int, subsampling_y : int, 
                             subsampling_z : int, block_id_x : int, block_id_y : int, block_id_z : int)
where reads(A_c, S, k), writes(S)
do

  var k_bounds = k.bounds.lo
  var rank = k[k_bounds]
  
  --column partition of A_c
  var A_c_bounds = A_c.bounds
  var n = A_c_bounds.hi.y - A_c_bounds.lo.y + 1
  var A_c_col_colors = ispace(int2d, {x = 1, y = n}) 
  var A_c_col_part = partition(equal, A_c, A_c_col_colors)  
  
  --column partition of S
  var S_bounds = S.bounds
  var S_col_part = partition(equal, S, A_c_col_colors)

  var y_points : int = ceil([double](ny)/[double](subsampling_y))
  var z_points : int = ceil([double](nz)/[double](subsampling_z))

  var block_index = block_id_z + 2*block_id_y + 2*2*block_id_x

  for col = 0,rank do

    var S_col = S_col_part[{x = 0, y = col}]
    var A_c_col = A_c_col_part[{x = 0, y = col}]

    var total_bounds = S_col.bounds
    var sub_reg_bounds = A_c_col.bounds
    var sub_reg_y = sub_reg_bounds.lo.y
    var sub_reg_x = sub_reg_bounds.lo.x

    var first_lower_sub_x_point = nx*block_id_x
    var first_lower_sub_y_point = ny*block_id_y
    var first_lower_sub_z_point = nz*block_id_z
  
    var first_upper_sub_x_point = nx*(block_id_x + 1) - 1
    var first_upper_sub_y_point = ny*(block_id_y + 1) - 1
    var first_upper_sub_z_point = nz*(block_id_z + 1) - 1
  
    var lower_sub_x = first_lower_sub_x_point
    var upper_sub_x = lower_sub_x + subsampling_x
    var lower_sub_y = first_lower_sub_y_point
    var upper_sub_y = lower_sub_y + subsampling_y
    var lower_sub_z = first_lower_sub_z_point
    var upper_sub_z = lower_sub_z + subsampling_z

    for i = 0, nx do

      var x_coord : int = i + block_id_x*nx 
      var x : double = dx*x_coord

      --check to see if subsampling interval has changed
      if x_coord >= upper_sub_x then
        lower_sub_x = lower_sub_x + subsampling_x
        upper_sub_x = upper_sub_x + subsampling_x
      end

      --check for points on the end of the boundary
      if x_coord == ((block_id_x + 1)*nx - 1) then
        upper_sub_x = x_coord
        lower_sub_x = upper_sub_x - subsampling_x
      end

      for j = 0, ny do

        var y_coord : int = j + block_id_y*ny
        var y : double = dy*y_coord

        if y_coord >= upper_sub_y then
          lower_sub_y = lower_sub_y + subsampling_y
          upper_sub_y = upper_sub_y + subsampling_y
        end

        --check for points on the end of the boundary
        if y_coord == ((block_id_y + 1)*ny - 1) then
          upper_sub_y = y_coord
          lower_sub_y = upper_sub_y - subsampling_y
        end

        for k = 0, nz do

          var z_coord : int = k + block_id_z*nz
          var z : double = dz*z_coord

        if z_coord >= upper_sub_z then
          lower_sub_z = lower_sub_z + subsampling_z
          upper_sub_z = upper_sub_z + subsampling_z
        end

          --check for points on the end of the boundary
          if z_coord == ((block_id_z + 1)*nz - 1) then
            upper_sub_z = z_coord
            lower_sub_z = upper_sub_z - subsampling_z
          end

          var x_0 : double = dx*lower_sub_x
          var y_0 : double = dy*lower_sub_y
          var z_0 : double = dz*lower_sub_z
          var x_1 : double = dx*upper_sub_x
          var y_1 : double = dy*upper_sub_y
          var z_1 : double = dz*upper_sub_z
          var x_01 : double = x_0 - x_1
          var y_01 : double = y_0 - y_1
          var z_01 : double = z_0 - z_1
          var x_10 : double = -1.0*x_01
          var y_10 : double = -1.0*y_01
          var z_10 : double = -1.0*z_01

          --adjust coords to be local
          var lower_sub_x_adj = lower_sub_x - nx*block_id_x
          var upper_sub_x_adj = upper_sub_x - nx*block_id_x
          var lower_sub_y_adj = lower_sub_y - ny*block_id_y
          var upper_sub_y_adj = upper_sub_y - nx*block_id_y
          var lower_sub_z_adj = lower_sub_z - nz*block_id_z
          var upper_sub_z_adj = upper_sub_z - nz*block_id_z

          --find neighboring points for subsampled cube
          var c_000_coord : int = lower_sub_z_adj/(subsampling_z) + z_points*(lower_sub_y_adj/(subsampling_y)) + z_points*y_points*(lower_sub_x_adj/(subsampling_x))
          var c_000 = A_c_col[{x = sub_reg_x + c_000_coord, y = sub_reg_y}] 

          var c_001_coord : int = upper_sub_z_adj/(subsampling_z) + z_points*(lower_sub_y_adj/(subsampling_y)) + z_points*y_points*(lower_sub_x_adj/(subsampling_x))
          var c_001 = A_c_col[{x = sub_reg_x + c_001_coord, y = sub_reg_y}]

          var c_010_coord : int = lower_sub_z_adj/subsampling_z + z_points*(upper_sub_y_adj/subsampling_y) + z_points*y_points*(lower_sub_x_adj/subsampling_x)
          var c_010 = A_c_col[{x = sub_reg_x + c_010_coord, y = sub_reg_y}]

          var c_011_coord : int = upper_sub_z_adj/subsampling_z + z_points*(upper_sub_y_adj/subsampling_y) + z_points*y_points*(lower_sub_x_adj/subsampling_x)
          var c_011 = A_c_col[{x = sub_reg_x + c_011_coord, y = sub_reg_y}]

          var c_100_coord : int = lower_sub_z_adj/subsampling_z + z_points*(lower_sub_y_adj/subsampling_y) + z_points*y_points*(upper_sub_x_adj/subsampling_x)
          var c_100 = A_c_col[{x = sub_reg_x + c_100_coord, y = sub_reg_y}]

          var c_101_coord : int = upper_sub_z_adj/subsampling_z + z_points*(lower_sub_y_adj/subsampling_y) + z_points*y_points*(upper_sub_x_adj/subsampling_x)
          var c_101 = A_c_col[{x = sub_reg_x + c_101_coord, y = sub_reg_y}]

          var c_110_coord : int = lower_sub_z_adj/subsampling_z + z_points*(upper_sub_y_adj/subsampling_y) + z_points*y_points*(upper_sub_x_adj/subsampling_x)
          var c_110 = A_c_col[{x = sub_reg_x + c_110_coord, y = sub_reg_y}]

          var c_111_coord : int = upper_sub_z_adj/subsampling_z + z_points*(upper_sub_y_adj/subsampling_y) + z_points*y_points*(upper_sub_x_adj/subsampling_x)
          var c_111 = A_c_col[{x = sub_reg_x + c_111_coord, y = sub_reg_y}]

          --compute coefficients for interpolated function
          var a_0 = (c_000*x_1*y_1*z_1)/(x_01*y_01*z_10) + (c_001*x_1*y_1*z_0)/(x_01*y_01*z_01) + (c_010*x_1*y_0*z_1)/(x_01*y_01*z_01)
                    + (c_011*x_1*y_0*z_0)/(x_01*y_01*z_10) + (c_100*x_0*y_1*z_1)/(x_01*y_01*z_01) + (c_101*x_0*y_1*z_0)/(x_01*y_01*z_10)
                    + (c_110*x_0*y_0*z_1)/(x_01*y_01*z_10) + (c_111*x_0*y_0*z_0)/(x_01*y_01*z_01)

          var a_1 = (c_000*y_1*z_1)/(x_01*y_01*z_01) + (c_001*y_1*z_0)/(x_01*y_01*z_10) + (c_010*y_0*z_1)/(x_01*y_01*z_10)
                    + (c_011*y_0*z_0)/(x_01*y_01*z_01) + (c_100*y_1*z_1)/(x_01*y_01*z_10) + (c_101*y_1*z_0)/(x_01*y_01*z_01)
                    + (c_110*y_0*z_1)/(x_01*y_01*z_01) + (c_111*y_0*z_0)/(x_01*y_01*z_10)

          var a_2 = (c_000*x_1*z_1)/(x_01*y_01*z_01) + (c_001*x_1*z_0)/(x_01*y_01*z_10) + (c_010*x_1*z_1)/(x_01*y_01*z_10)
                    + (c_011*x_1*z_0)/(x_01*y_01*z_01) + (c_100*x_0*z_1)/(x_01*y_01*z_10) + (c_101*x_0*z_0)/(x_01*y_01*z_01)
                    + (c_110*x_0*z_1)/(x_01*y_01*z_01) + (c_111*x_0*z_0)/(x_01*y_01*z_10)

          var a_3 = (c_000*x_1*y_1)/(x_01*y_01*z_01) + (c_001*x_1*y_1)/(x_01*y_01*z_10) + (c_010*x_1*y_0)/(x_01*y_01*z_10)
                    + (c_011*x_1*y_0)/(x_01*y_01*z_01) + (c_100*x_0*y_1)/(x_01*y_01*z_10) + (c_101*x_0*y_1)/(x_01*y_01*z_01)
                    + (c_110*x_0*y_0)/(x_01*y_01*z_01) + (c_111*x_0*y_0)/(x_01*y_01*z_10)

          var a_4 = (c_000*z_1)/(x_01*y_01*z_10) + (c_001*z_0)/(x_01*y_01*z_01) + (c_010*z_1)/(x_01*y_01*z_01)
                    + (c_011*z_0)/(x_01*y_01*z_10) + (c_100*z_1)/(x_01*y_01*z_01) + (c_101*z_0)/(x_01*y_01*z_10)
                    + (c_110*z_1)/(x_01*y_01*z_10) + (c_111*z_0)/(x_01*y_01*z_01)

          var a_5 = (c_000*y_1)/(x_01*y_01*z_10) + (c_001*y_1)/(x_01*y_01*z_01) + (c_010*y_0)/(x_01*y_01*z_01)
                    + (c_011*y_0)/(x_01*y_01*z_10) + (c_100*y_1)/(x_01*y_01*z_01) + (c_101*y_1)/(x_01*y_01*z_10)
                    + (c_110*y_0)/(x_01*y_01*z_10) + (c_111*y_0)/(x_01*y_01*z_01)

          var a_6 = (c_000*x_1)/(x_01*y_01*z_10) + (c_001*x_1)/(x_01*y_01*z_01) + (c_010*x_1)/(x_01*y_01*z_01)
                    + (c_011*x_1)/(x_01*y_01*z_10) + (c_100*x_0)/(x_01*y_01*z_01) + (c_101*x_0)/(x_01*y_01*z_10)
                    + (c_110*x_0)/(x_01*y_01*z_10) + (c_111*x_0)/(x_01*y_01*z_01)

          var a_7 = (c_000)/(x_01*y_01*z_01) + (c_001)/(x_01*y_01*z_10) + (c_010)/(x_01*y_01*z_10)
                    + (c_011)/(x_01*y_01*z_01) + (c_100)/(x_01*y_01*z_10) + (c_101)/(x_01*y_01*z_01)
                    + (c_110)/(x_01*y_01*z_01) + (c_111)/(x_01*y_01*z_10)

          --construct interpolated function
          S_col[{x = total_bounds.lo.x + (k + j*nz + i*nz*ny), y = total_bounds.lo.y}] = a_0 + a_1*x + a_2*y + a_3*z 
                                                                                 + a_4*x*y + a_5*x*z + a_6*y*z + a_7*x*y*z
       end --end z-dir loop
     
       --reset z_coord values 
       lower_sub_z = first_lower_sub_z_point
       upper_sub_z = lower_sub_z + subsampling_z

      end  --end y-dir loop

      --reset y_coord values 
      lower_sub_y = first_lower_sub_y_point
      upper_sub_y = lower_sub_y + subsampling_y

    end  --end x-dir loop

    --reset x_coord values 
    lower_sub_x = first_lower_sub_x_point
    upper_sub_x = lower_sub_x + subsampling_x

  end  --end rank loop

end

--task to reconstruct the 3D domain solution at a specific time step
task id_exp.reconstruct_solution_time(A_c : region(ispace(int2d), double), P_col : region(ispace(int2d), double), domain : region(ispace(int3d), field_space), k : region(ispace(int2d), int))
where reads(A_c, P_col, k), writes(domain.[qoi_fld])
do

  var k_bounds = k.bounds.lo
  var rank = k[k_bounds]
  
  --column partition of A_c
  var A_c_bounds = A_c.bounds
  var A_c_row_dim = {x = A_c_bounds.hi.x - A_c_bounds.lo.x + 1, y = 1}
  var A_c_row_colors = ispace(int2d, A_c_row_dim)
  var A_c_row_part = partition(equal, A_c, A_c_row_colors)

  --find the 1D domain reconstruction of the solution
  var domain_1d = region(ispace(int1d, A_c_bounds.hi.x - A_c_bounds.lo.x + 1), double)
  fill(domain_1d, 0.0)

  for k in domain_1d do
  var A_c_row = A_c_row_part[{x = k, y = 0}]
    for j = P_col.bounds.lo.x, (P_col.bounds.lo.x + rank) do
      domain_1d[k] = domain_1d[k] + P_col[{x = j, y = P_col.bounds.lo.y}]*A_c_row[{x = A_c_row.bounds.lo.x, y = A_c_row.bounds.lo.y + j - P_col.bounds.lo.x}] 
    end
  end

  --go from 1d domain to 3d domain
  var domain_bounds = domain.bounds
  var nz = domain_bounds.hi.z - domain_bounds.lo.z + 1
  var ny = domain_bounds.hi.y - domain_bounds.lo.y + 1
  
  for i in domain do
    domain[i].[qoi_fld] = domain_1d[(i.z - domain_bounds.lo.z) + (i.y - domain_bounds.lo.y)*nz + (i.x - domain_bounds.lo.x)*nz*ny]
  end

end

--task to copy contents of one region to another
task id_exp.copy_qoi(source : region(ispace(int2d), double), target : region(ispace(int2d), double))
where reads(source, target), writes(target)
do
  copy(source, target)
end

--task to find the Frobenius norm of the difference between exact and approximate solutions
task id_exp.error_frob_norm_3d(domain_restart : region(ispace(int3d), field_space), domain_no_restart : region(ispace(int3d), field_space))
where reads(domain_restart.[qoi_fld], domain_no_restart.[qoi_fld])
do

  var rs_bounds = domain_restart.bounds
  var x_dim = rs_bounds.hi.x - rs_bounds.lo.x + 1
  var y_dim = rs_bounds.hi.y - rs_bounds.lo.y + 1
  var z_dim = rs_bounds.hi.z - rs_bounds.lo.z + 1

  var no_rs_bounds = domain_no_restart.bounds

  --find the difference between the two solutions
  var domain_diff = region(ispace(int3d, {x = x_dim, y = y_dim, z = z_dim}), double)

  for i in domain_diff do
    domain_diff[i] = domain_restart[rs_bounds.lo + i].[qoi_fld] - domain_no_restart[no_rs_bounds.lo + i].[qoi_fld] 
  end

  --compute the Frobenius norm of the difference matrix
  var norm_frob : double = 0.0
 
  for i in domain_diff do
    norm_frob = norm_frob + domain_diff[i]*domain_diff[i]
  end

  --compute the Frobenius norm of the "exact" solution (used for normalization)
  var norm_frob_A : double = 0.0
 
  for i in domain_no_restart do
    norm_frob_A = norm_frob_A + domain_no_restart[i].[qoi_fld]*domain_no_restart[i].[qoi_fld]
  end

  return sqrt(norm_frob)/sqrt(norm_frob_A)

end



return id_exp

end
