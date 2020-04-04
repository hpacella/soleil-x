import "regent"

local c = regentlib.c
local std = terralib.includec("stdlib.h")
local sqrt = regentlib.sqrt(double)
local abs = regentlib.fabs(double)
--local blas = require("BLAS_functions")
local DBL_EPSILON = 2.220446e-16
local fmod = regentlib.fmod(double)

local id_la = {}

--task to generate a random i.i.d. gaussian matrix
task id_la.make_gaussian_matrix(omega : region(ispace(int2d), double))
where reads(omega), writes(omega) 
do
  for i in omega do
    omega[i] = [double](std.drand48())
  end

end

--task to transpose a matrix region
__demand(__inline)
task transpose_region(original : region(ispace(int2d), double), transpose : region(ispace(int2d), double))
where reads(original), writes(transpose)
do

  var orig_bounds = original.bounds
  var m_orig = orig_bounds.hi.y - orig_bounds.lo.y + 1
  var n_orig = orig_bounds.hi.x - orig_bounds.lo.x + 1

  var trans_bounds = transpose.bounds
  var m_trans = trans_bounds.hi.y - trans_bounds.lo.y + 1
  var n_trans = trans_bounds.hi.x - trans_bounds.lo.x + 1
  
  var x_shift = orig_bounds.lo.x
  var y_shift = orig_bounds.lo.y

  if m_orig ~= n_trans or n_orig ~= m_trans then
    c.printf("Dimensions do not match!\n")
    c.abort()
  end

  for i in original do
    transpose[ {x = i.y - y_shift, y = i.x + y_shift}] = original[i]
  end

end

--task to (left or right) multiply a matrix and a diagonal matrix
__demand(__inline)
task diag_matrix_mult(diag_matrix : region(ispace(int1d), double), matrix : region(ispace(int2d), double), dir :int)
where reads(matrix, diag_matrix), writes(matrix)
do

  var m_bounds = matrix.bounds
  var s_bounds = diag_matrix.bounds

  if dir == 0 then --left multiplication by diagonal matrix 
    for i in matrix do
      matrix[i] = 1/diag_matrix[i.x - m_bounds.lo.x + s_bounds.lo]*matrix[i]
    end
      
  elseif dir == 1 then --right multiplication by diagonal matrix
    for i in matrix do
      matrix[i] = 1/diag_matrix[i.y - m_bounds.lo.y + s_bounds.lo]*matrix[i]
    end
      
  else
    c.printf("Provide valid direction.\n")
  end
      
end
      
--task to find pivot column of a matrix
__demand(__inline)
task find_column_pivot(domain : region(ispace(int2d), double))
where reads(domain)
do

  var sq_entry : double = 0.0

  for i in domain do
    sq_entry = sq_entry + domain[i]*domain[i]
  end
  
  return sq_entry

end

--task to compute the vector 2-norm
__demand(__inline)
task vector_two_norm(vector : region(ispace(int2d), double))
where reads(vector)
do

  var two_norm : double = 0.0

  for i in vector do
    two_norm += vector[i]*vector[i]
  end

  return sqrt(two_norm)

end

__demand(__inline)
task swap_cols(column_1 : region(ispace(int2d), double), column_2 : region(ispace(int2d), double), m : int, imax : int)
where
reads(column_1, column_2), writes(column_1, column_2)
do
  
  var temp_array_1 = region(ispace(int1d, m), double)
  var temp_array_2 = region(ispace(int1d, m), double)

  --copy overwritten column to temporary array
  for i in column_1 do
    temp_array_1[i.x] = column_1[i]
  end
  
  for i in column_2 do
    temp_array_2[i.x] = column_2[i]
  end

  --overwrite original column with pivot column
  for i in column_1 do
    column_1[i] = temp_array_2[i.x]
  end
  
  --copy temp to new column
  for i in column_2 do
    column_2[i] = temp_array_1[i.x] 
  end
  
end

__demand(__inline)
task R_entry_update(column_k : region(ispace(int2d), double), R_kk : double, column_i : region(ispace(int2d), double), i : int)
where reads(column_k, column_i)
do

  var R_entry : double = 0.0
  var column_i_ind = column_i.bounds.lo.y

  for j in column_k do
    R_entry = R_entry + (1.0/R_kk*column_k[j])*column_i[{x = j.x, y = column_i_ind}] 
  end

  return R_entry

end 

__demand(__inline)
task column_update(column_j : region(ispace(int2d), double), R_kj : double, R_kk : double, column_k : region(ispace(int2d), double))
where reads(column_k), reduces -(column_j) 
do

  var column_k_ind = column_k.bounds.lo.y

  for i in column_j do
    column_j[i] -= R_kj*(1.0/R_kk*column_k[{x = i.x, y = column_k_ind}])
  end

end

--task to perform DGEMM
__demand(__inline)
task dgemm(C_m : int, C_n : int, C_k : int, C : region(ispace(int2d), double), A : region(ispace(int2d), double), B : region(ispace(int2d), double))
where reads(A, B, C), writes(C)
do

  for k = 0, C_k do
    for j = 0, C_n do
      for i = 0, C_m do
        C[{x = C.bounds.lo.x + i, y = C.bounds.lo.y + j}] = C[{x = C.bounds.lo.x + i, y = C.bounds.lo.y + j}] 
              + A[{x = A.bounds.lo.x + i, y = A.bounds.lo.y + k}]*B[{x = B.bounds.lo.x + k, y = B.bounds.lo.y + j}]
      end
    end
  end

end

--task to perform backwards substitution for a system with an upper triangular matrix
__demand(__inline)
task upper_tri_back_sub(C_m : int, C_n : int, C_k : int, C : region(ispace(int2d), double), U : region(ispace(int2d), double), B : region(ispace(int2d), double))
where reads(C, U, B), writes(C)
do

  var inner_prod : double
  
  --want to solve the linear system C = U\B
  --C_n = the number of columns of B
  --C_m = the number of rows of C

  for i = 0, C_n do --traverse columns of B
    
    C[{x = C.bounds.hi.x, y = C.bounds.lo.y + i}] = B[{x = B.bounds.hi.x, y = B.bounds.lo.y + i}]/U[{x = U.bounds.hi.x, y = U.bounds.hi.y}]

    for j = (C_m - 2), -1, -1  do  --traverse rows of C
      
      for k = (j + 1), C_k do

        inner_prod += U[{x = U.bounds.lo.x + j, y = U.bounds.lo.y + k}]*C[{x = C.bounds.lo.x + k, y = C.bounds.lo.y + i}]

      end

      C[{x = C.bounds.lo.x + j, y = C.bounds.lo.y + i}] = (B[{x = B.bounds.lo.y + j, y = B.bounds.lo.y + i}] - inner_prod)/U[{x = U.bounds.lo.x + j, y = U.bounds.lo.y + j}]
      
      inner_prod = 0.0

    end
  end

end

--task to find pivoted QR using modified Gram-Schmidt
task id_la.complete_ID(domain : region(ispace(int2d), double), global_pivs : region(ispace(int2d),int), tol : region(ispace(int2d), double), 
                       k_rank : region(ispace(int2d), int), P_final : region(ispace(int2d), double))
where reads(domain, global_pivs, tol, P_final), writes(domain, global_pivs, k_rank, P_final)
do

  --matrix dimensions
  var bounds = domain.bounds
  var global_pivs_x = global_pivs.bounds.lo.x
  var m = bounds.hi.x - bounds.lo.x + 1  --no. of rows
  var n = bounds.hi.y - bounds.lo.y + 1  --no. of columns
  var maxrnk = min(m,n)
  var column_sums = region(ispace(int1d, n), double)
  var tol_idx = tol.bounds.lo
  var k_bounds = k_rank.bounds.lo 

  --copy domain into a temp region to be operated on
  var domain_copy = region(ispace(int2d, {x = m, y = n}), double)
  for i in domain do
    domain_copy[{x = i.x - bounds.lo.x, y = i.y - bounds.lo.y}] = domain[i]
  end

  --initialize pivot vector
  var pivot_cols = region(ispace(int2d, {x = 1, y = n}), int)
  for i in pivot_cols do
    pivot_cols[i] = i.y
  end

  --create column partition of the domain and R
  var R = region(ispace(int2d, {x = n, y = n}), double)
  var R_bounds = R.bounds
  var column_space = ispace(int2d, {x = 1, y = n}) 
  var column_part = partition(equal, domain, column_space)
  var copy_column_part = partition(equal, domain_copy, column_space)
  var R_column_part = partition(equal, R, column_space)

  --find inital pivot column
  var cmax : double = 0.0
  var imax : int = 0
 
  for i in column_space do

    column_sums[i.y] = find_column_pivot(copy_column_part[i])
    
    if column_sums[i.y] > cmax then
      cmax = column_sums[i.y]
      imax = i.y
    end
  end

  var rank : int = 0
    
  for k = 0, maxrnk do

    --perform pivot
    var tmp_pivot = pivot_cols[{x = 0, y = k}]
    pivot_cols[{x = 0, y = k}] = pivot_cols[{x = 0, y = imax}]
    pivot_cols[{x = 0, y = imax}] = tmp_pivot
    
    var tmp_global = global_pivs[{x = global_pivs_x, y = k}]
    global_pivs[{x = global_pivs_x, y = k}] = global_pivs[{x = global_pivs_x, y = imax}]
    global_pivs[{x = global_pivs_x, y = imax}] = tmp_global

    swap_cols(copy_column_part[{x = 0, y = k}], copy_column_part[{x = 0, y = imax}], m, imax)
    swap_cols(column_part[{x = 0, y = k}], column_part[{x = 0, y = imax}], m, imax)
    swap_cols(R_column_part[{x = 0, y = k}], R_column_part[{x = 0, y = imax}], m, imax)
    
    var tmp_c = column_sums[k]
    column_sums[k] = column_sums[imax]
    column_sums[imax] = tmp_c    
    
    --Gram-Schmidt sweep
    var R_kk : double = vector_two_norm(copy_column_part[{x = 0, y = k}])
    R[{x = R_bounds.lo.x + k, y = R_bounds.lo.y + k}] = R_kk
   
    for i = (k+1), n do
      R[{x = R_bounds.lo.x + k, y = R_bounds.lo.y + i}] = R_entry_update(copy_column_part[{x = 0, y = k}], R_kk, 
                                                          copy_column_part[{x = 0, y = i}], i)
    end

    for j = (k+1), n do
      column_update(copy_column_part[{x = 0, y = j}], R[{x = R_bounds.lo.x + k, y = R_bounds.lo.y + j}], 
                    R_kk, copy_column_part[{x = 0, y = k}]) 
    end

    --determine next pivot
    if k < (maxrnk - 1) then
      cmax = 0.0
      for i = (k+1), n do
        column_sums[i] = find_column_pivot(copy_column_part[{x = 0, y = i}])
        if column_sums[i] > cmax then
          cmax = column_sums[i]
          imax = i
        end
      end
    end

    --check for convergence
    if abs(cmax) < tol[tol_idx] then
      rank = k + 1
      break
    end

  if k == maxrnk - 1 then
    rank = k + 1
  end

  end  --end for loop
  
  k_rank[k_bounds] = rank
 
  --compute P

  --find R_11 and R_12 partitions
  var R_11 = region(ispace(int2d, {x = rank, y = rank}), double)
  var R_12 = region(ispace(int2d, {x = rank, y = (n - rank)}), double)

  --extract R_11
  var R_11_bounds = R_11.bounds

  for i = R_11_bounds.lo.x, (R_11_bounds.lo.x + rank) do
    for j = R_11_bounds.lo.y, (R_11_bounds.lo.y + rank) do
    R_11[ {x = i, y = j}] = R[{x = i - R_11_bounds.lo.x + R_bounds.lo.x, 
                            y = j - R_11_bounds.lo.y + R_bounds.lo.y}]
    end
  end
  
  --extract R_12
  var R_12_bounds = R_12.bounds

  for i = (R_12_bounds.lo.x), (R_12_bounds.hi.x + 1) do
    for j = (R_12_bounds.lo.y), (R_12_bounds.hi.y + 1) do
      R_12[ {x = i, y = j}] = R[{x = i - R_12_bounds.lo.x + R_bounds.lo.x, 
                            y = j - R_12_bounds.lo.y + R_bounds.lo.y + rank}]
    end
  end

  --construct P(:,I)
  var P = region(ispace(int2d, {x = rank, y = n}), double) 
  var P_bounds = P.bounds
  fill(P, 0.0)

  for i = P_bounds.lo.x, P_bounds.lo.x + rank do
    for j = P_bounds.lo.y, P_bounds.hi.y + 1 do 
      if (i - P_bounds.lo.x) == (j - P_bounds.lo.y) then
        P[{x = i, y = j}] = 1.0
      end
    end
  end

  --construct permutation matrix Perm
  var Perm = region(ispace(int2d, {x = n, y = n}), double)
  fill(Perm, 0.0) 
  for i in pivot_cols do
    Perm[{x = Perm.bounds.lo.x + pivot_cols[i], y = Perm.bounds.lo.y + i.y}] = 1.0
  end

  if rank ~= n then
 
    --T = R_11^(-1)*R_12 (or identity matrix if R_11 == R)
    var T = region(ispace(int2d, {x = rank, y = (n - rank)}), double)
    upper_tri_back_sub(rank, (n - rank), rank, T, R_11, R_12)

    --copy T into P
    var P_bounds = P.bounds
    var T_bounds = T.bounds

    for i = (P_bounds.lo.x), (P_bounds.lo.x + rank) do
      for j = (P_bounds.hi.y - (n - rank) + 1), (P_bounds.hi.y + 1) do
      P[{x = i, y = j}] = T[{x = i - P_bounds.lo.x, y = j - (P_bounds.hi.y - (n - rank) + 1)}]
      end
    end
    
    --transpose the Perm matrix
    var Perm_T = region(ispace(int2d, {x = n, y = n}), double)
    transpose_region(Perm, Perm_T)

    --compute final coefficient matrix    
    dgemm(rank, n, n, P_final, P, Perm_T)

  else --edge case when R is full rank, P(:,I) is the identity matrix

    --transpose the Perm matrix
    var Perm_T = region(ispace(int2d, {x = n, y = n}), double)
    transpose_region(Perm, Perm_T)

    --compute final coefficient matrix
    dgemm(rank, n, n, P_final, P, Perm_T)

  end


end

--task to perform MGSQR pivot of the columns of the matrix T_i
task id_la.complete_ID_on_T_i(T_i : region(ispace(int2d), double), T_i_pivot_cols : region(ispace(int2d),int), time_step : int,
                     target_rank : int, P_final : region(ispace(int2d), double), global_piv : region(ispace(int2d), int),
                      A_c : region(ispace(int2d), double), n_total : int, no_tstep_intervals : int)
where reads(T_i, T_i_pivot_cols, P_final), writes(T_i, T_i_pivot_cols, P_final, global_piv, A_c)
do

  --matrix dimensions
  var bounds = T_i.bounds
  var pivot_bounds = T_i_pivot_cols.bounds
  var pivot_cols_x = pivot_bounds.lo.x
  var m = bounds.hi.x - bounds.lo.x + 1  --no. of rows
  var n = bounds.hi.y - bounds.lo.y + 1  --no. of columns (=no. of time steps in each subinterval)
  var column_sums = region(ispace(int1d, n), double)

  --create R region and partition into columns
  var R = region(ispace(int2d, {x = target_rank, y = n}), double)
  var column_space = ispace(int2d, {x = 1, y = n}) 
  var R_column_part = partition(equal, R, column_space)

  --copy T_i into a temp region to be operated on
  var T_i_copy = region(ispace(int2d, {x = m, y = n}), double)
  for i in T_i do
    T_i_copy[{x = i.x - bounds.lo.x, y = i.y - bounds.lo.y}] = T_i[i]
  end

  --create column partition of T_i
  var column_part = partition(equal, T_i, column_space)
  var copy_column_part = partition(equal, T_i_copy, column_space)

  --initialize the pivot cols
  for i in T_i_pivot_cols do
    T_i_pivot_cols[i] = i.y
  end

  --find inital pivot column
  var cmax : double = 0.0
  var imax : int = 0
 
  for i in column_space do

    column_sums[i.y] = find_column_pivot(copy_column_part[i])
    
    if column_sums[i.y] > cmax then
      cmax = column_sums[i.y]
      imax = i.y
    end
  end

  for k = 0,target_rank do

    --perform pivot on T_i, T_i_copy
    var tmp_pivot = T_i_pivot_cols[{x = pivot_cols_x, y = k}]
    T_i_pivot_cols[{x = pivot_cols_x, y = k}] = T_i_pivot_cols[{x = pivot_cols_x, y = imax}]
    T_i_pivot_cols[{x = pivot_cols_x, y = imax}] = tmp_pivot
    
    swap_cols(copy_column_part[{x = 0, y = k}], copy_column_part[{x = 0, y = imax}], m, imax)
    swap_cols(column_part[{x = 0, y = k}], column_part[{x = 0, y = imax}], m, imax)
    swap_cols(R_column_part[{x = 0, y = k}], R_column_part[{x = 0, y = imax}], m, imax)
    
    var tmp_c = column_sums[k]
    column_sums[k] = column_sums[imax]
    column_sums[imax] = tmp_c    
    
    --Gram-Schmidt sweep
    var R_kk : double = vector_two_norm(copy_column_part[{x = 0, y = k}])
    R[{x = k, y = k}] = R_kk

    for j = (k+1),n do
      var R_kj = R_entry_update(copy_column_part[{x = 0, y = k}], R_kk, copy_column_part[{x = 0, y = j}], j)
      R[{x = k, y = j}] = R_kj
      column_update(copy_column_part[{x = 0, y = j}], R_kj, R_kk, copy_column_part[{x = 0, y = k}]) 
    end

    --determine next pivot
    if k < (target_rank - 1) then
      cmax = 0.0
      for i = (k+1),n do
        column_sums[i] = find_column_pivot(copy_column_part[{x = 0, y = i}])
        if column_sums[i] > cmax then
          cmax = column_sums[i]
          imax = i 
        end
      end
    
    end

  end  --end for loop

  --convert local pivots to global indexing
  var offset : int = target_rank*((time_step + 1)/(n_total/no_tstep_intervals) - 1)
  for i = 0, target_rank do
    global_piv[{x = -1*T_i_pivot_cols.bounds.lo.x + global_piv.bounds.lo.x, y = i - T_i_pivot_cols.bounds.lo.y + 
                 global_piv.bounds.lo.y + offset}] = T_i_pivot_cols[{x = T_i_pivot_cols.bounds.lo.x, 
                 y = T_i_pivot_cols.bounds.lo.y + i}] + (offset/target_rank)*(n_total/no_tstep_intervals)
  end
  
  --compute P
 
  --find R_11 and R_12 partitions
  var R_11 = region(ispace(int2d, {x = target_rank, y = target_rank}), double)
  var R_12 = region(ispace(int2d, {x = target_rank, y = (n - target_rank)}), double)

  --extract R_11
  var R_11_bounds = R_11.bounds

  for i = R_11_bounds.lo.x, (R_11_bounds.lo.x + target_rank) do
    for j = R_11_bounds.lo.y, (R_11_bounds.lo.y + target_rank) do
    R_11[ {x = i, y = j}] = R[{x = i - R_11_bounds.lo.x, 
                            y = j - R_11_bounds.lo.y}]
    end
  end
  
  --extract R_12
  var R_12_bounds = R_12.bounds

  for i = (R_12_bounds.lo.x), (R_12_bounds.hi.x + 1) do
    for j = (R_12_bounds.lo.y), (R_12_bounds.hi.y + 1) do
      R_12[ {x = i, y = j}] = R[{x = i - R_12_bounds.lo.x, 
                            y = j - R_12_bounds.lo.y + target_rank}]
    end
  end

  --construct P(:,I)
  var P = region(ispace(int2d, {x = target_rank, y = n}), double) 
  var P_bounds = P.bounds
  fill(P, 0.0)

  for i = P_bounds.lo.x, P_bounds.lo.x + target_rank do
    for j = P_bounds.lo.y, P_bounds.lo.y + target_rank do 
      if (i - P_bounds.lo.x) == (j - P_bounds.lo.y) then
        P[{x = i, y = j}] = 1.0
      end
    end
  end

  --construct permutation matrix Perm
  var Perm = region(ispace(int2d, {x = n, y = n}), double)
  fill(Perm, 0.0) 
  for i in T_i_pivot_cols do
    Perm[{x = Perm.bounds.lo.x + T_i_pivot_cols[i], y = Perm.bounds.lo.y + i.y}] = 1.0
  end

  if target_rank ~= n then
 
    --T = R_11^(-1)*R_12 (or identity matrix if R_11 == R)
    var T = region(ispace(int2d, {x = target_rank, y = (n - target_rank)}), double)
    upper_tri_back_sub(target_rank, (n - target_rank), target_rank, T, R_11, R_12)

    --copy T into P
    var P_bounds = P.bounds
    var T_bounds = T.bounds

    for i = (P_bounds.lo.x), (P_bounds.lo.x + target_rank) do
      for j = (P_bounds.hi.y - (n - target_rank) + 1), (P_bounds.hi.y + 1) do
      P[{x = i, y = j}] = T[{x = i - P_bounds.lo.x, y = j - (P_bounds.hi.y - (n - target_rank) + 1)}]
      end
    end

    --transpose the Perm matrix
    var Perm_T = region(ispace(int2d, {x = n, y = n}), double)
    transpose_region(Perm, Perm_T)

    --compute final coefficient matrix    
    dgemm(target_rank, n, n, P_final, P, Perm_T)

  else --edge case when R is full rank, P(:,I) is the identity matrix

    --transpose the Perm matrix
    var Perm_T = region(ispace(int2d, {x = n, y = n}), double)
    transpose_region(Perm, Perm_T)

    --compute final coefficient matrix
    dgemm(target_rank, n, n, P_final, P, Perm_T)

  end

  --copy local results to global A_c matrix
  var A_c_bounds = A_c.bounds
  var A_c_offset : int = target_rank*((time_step + 1)/(n_total/no_tstep_intervals) - 1)

  for j = bounds.lo.x, bounds.hi.x + 1 do
    for k = bounds.lo.y, bounds.lo.y + (target_rank + 1) do
      A_c[{x = j - bounds.lo.x + A_c_bounds.lo.x, y = k - bounds.lo.y + A_c_bounds.lo.y + A_c_offset}] = T_i[{x = j, y = k}]
    end
  end

end

--[[
--task to compute coefficient matrix P using BLAS subroutines
task id_la.compute_P_BLAS(R : region(ispace(int2d), double), k_rank : region(ispace(int2d), int), P_final : region(ispace(int2d), double),
                      pivot_columns : region(ispace(int2d), int), x_ind : int, y_ind : int)
where reads(R, k_rank, pivot_columns, P_final), writes(P_final)
do

  var bounds = R.bounds
  var n = bounds.hi.y - bounds.lo.y + 1

  --find R_11 and R_12 partitions
  var rank: int
  for i in k_rank do
    rank = k_rank[i]
  end
  c.printf("rank = %d\n", rank)
  c.printf("n = %d\n", n)
  var R_11 = region(ispace(int2d, {x = rank, y = rank}), double)
  var R_12 = region(ispace(int2d, {x = rank, y = (n - rank)}), double)

  --extract R_11
  var R_11_bounds = R_11.bounds

  for i = R_11_bounds.lo.x, (R_11_bounds.lo.x + rank) do
    for j = R_11_bounds.lo.y, (R_11_bounds.lo.y + rank) do
    R_11[ {x = i, y = j}] = R[{x = i - R_11_bounds.lo.x + bounds.lo.x, 
                            y = j - R_11_bounds.lo.y + bounds.lo.y}]
    end
  end
  
  --extract R_12
  var R_12_bounds = R_12.bounds

  for i = (R_12_bounds.lo.x), (R_12_bounds.hi.x + 1) do
    for j = (R_12_bounds.lo.y), (R_12_bounds.hi.y + 1) do
      R_12[ {x = i, y = j}] = R[{x = i - R_12_bounds.lo.x + bounds.lo.x, 
                            y = j - R_12_bounds.lo.y + bounds.lo.y + rank}]
    end
  end

  --construct P(:,I)
  var P = region(ispace(int2d, {x = rank, y = n}), double) 
  var P_bounds = P.bounds
  fill(P, 0.0)

  for i = P_bounds.lo.x, P_bounds.lo.x + rank do
    for j = P_bounds.lo.y, P_bounds.hi.y + 1 do 
      if (i - P_bounds.lo.x) == (j - P_bounds.lo.y) then
        P[{x = i, y = j}] = 1.0
      end
    end
  end

  --construct permutation matrix Perm
  var Perm = region(ispace(int2d, {x = n, y = n}), double)
  fill(Perm, 0.0) 
  for i in pivot_columns do
    Perm[{x = Perm.bounds.lo.x + pivot_columns[i], y = Perm.bounds.lo.y + i.y}] = 1.0
  end

  if rank ~= n then
 
    --T = R_11^(-1)*R_12 (or identity matrix if R_11 == R)
    var T = region(ispace(int2d, {x = rank, y = (n - rank)}), double)

    --compute the SVD of R_11
    var U = region(ispace(int2d, {x = rank, y = rank}), double)
    var U_T = region(ispace(int2d, {x = rank, y = rank}), double)
    var V_T = region(ispace(int2d, {x = rank, y = rank}), double)
    var S = region(ispace(int1d, rank), double)
    var R_pinv = region(ispace(int2d, {x = rank, y = rank}), double)

    var svd_work = region(ispace(int1d, 10*rank*rank), double) 
    var svd_iwork = region(ispace(int1d, 12*rank), int)

    blas.dgesvdx(R_11, U, V_T, S, rank, svd_work, svd_iwork, 0, 0) 

    --transpose U
    transpose_region(U, U_T)
  
    --multiply inv(S)*U^T
    diag_matrix_mult(S, U_T, 0)

    --multiply V*(inv(S)*U^T) = R_pinv 
    blas.dgemm(V_T, U_T, R_pinv, 1.0, 0.0, 1, 0, rank, rank, rank, 0, 0, 0, 0, 0, 0)
  
    --find T = R_pinv*R_12
    blas.dgemm(R_pinv, R_12, T, 1.0, 0.0, 0, 0, (n-rank), rank, rank, 0, 0, 0, 0, 0, 0)
  
    --copy T into P
    var P_bounds = P.bounds
    var T_bounds = T.bounds

    for i = (P_bounds.lo.x), (P_bounds.lo.x + rank) do
      for j = (P_bounds.hi.y - (n - rank) + 1), (P_bounds.hi.y + 1) do
      P[{x = i, y = j}] = T[{x = i - P_bounds.lo.x, y = j - (P_bounds.hi.y - (n - rank) + 1)}]
      end
    end
    
    --compute final coefficient matrix
    blas.dgemm(P, Perm, P_final, 1.0, 0.0, 0, 1, n, rank, n, 0, 0, 0, 0, x_ind, y_ind)

  else --edge case when R is full rank, P(:,I) is the identity matrix

    --compute final coefficient matrix
    blas.dgemm(P, Perm, P_final, 1.0, 0.0, 0, 1, n, rank, n, 0, 0, 0, 0, x_ind, y_ind)

  end

end
]]--

--task to compute coefficient matrix P
task id_la.compute_P(R : region(ispace(int2d), double), k_rank : region(ispace(int2d), int), P_final : region(ispace(int2d), double),
                      pivot_columns : region(ispace(int2d), int), x_ind : int, y_ind : int)
where reads(R, k_rank, pivot_columns, P_final), writes(P_final)
do

  var bounds = R.bounds
  var n = bounds.hi.y - bounds.lo.y + 1

  --find R_11 and R_12 partitions
  var rank: int
  for i in k_rank do
    rank = k_rank[i]
  end
  var R_11 = region(ispace(int2d, {x = rank, y = rank}), double)
  var R_12 = region(ispace(int2d, {x = rank, y = (n - rank)}), double)

  --extract R_11
  var R_11_bounds = R_11.bounds

  for i = R_11_bounds.lo.x, (R_11_bounds.lo.x + rank) do
    for j = R_11_bounds.lo.y, (R_11_bounds.lo.y + rank) do
    R_11[ {x = i, y = j}] = R[{x = i - R_11_bounds.lo.x + bounds.lo.x, 
                            y = j - R_11_bounds.lo.y + bounds.lo.y}]
    end
  end
  
  --extract R_12
  var R_12_bounds = R_12.bounds

  for i = (R_12_bounds.lo.x), (R_12_bounds.hi.x + 1) do
    for j = (R_12_bounds.lo.y), (R_12_bounds.hi.y + 1) do
      R_12[ {x = i, y = j}] = R[{x = i - R_12_bounds.lo.x + bounds.lo.x, 
                            y = j - R_12_bounds.lo.y + bounds.lo.y + rank}]
    end
  end

  --construct P(:,I)
  var P = region(ispace(int2d, {x = rank, y = n}), double) 
  var P_bounds = P.bounds
  fill(P, 0.0)

  for i = P_bounds.lo.x, P_bounds.lo.x + rank do
    for j = P_bounds.lo.y, P_bounds.hi.y + 1 do 
      if (i - P_bounds.lo.x) == (j - P_bounds.lo.y) then
        P[{x = i, y = j}] = 1.0
      end
    end
  end

  --construct permutation matrix Perm
  var Perm = region(ispace(int2d, {x = n, y = n}), double)
  fill(Perm, 0.0) 
  for i in pivot_columns do
    Perm[{x = Perm.bounds.lo.x + pivot_columns[i], y = Perm.bounds.lo.y + i.y}] = 1.0
  end

  if rank ~= n then
 
    --T = R_11^(-1)*R_12 (or identity matrix if R_11 == R)
    var T = region(ispace(int2d, {x = rank, y = (n - rank)}), double)
    upper_tri_back_sub(rank, (n - rank), rank, T, R_11, R_12)

    --copy T into P
    var P_bounds = P.bounds
    var T_bounds = T.bounds

    for i = (P_bounds.lo.x), (P_bounds.lo.x + rank) do
      for j = (P_bounds.hi.y - (n - rank) + 1), (P_bounds.hi.y + 1) do
      P[{x = i, y = j}] = T[{x = i - P_bounds.lo.x, y = j - (P_bounds.hi.y - (n - rank) + 1)}]
      end
    end
    
    --transpose the Perm matrix
    var Perm_T = region(ispace(int2d, {x = n, y = n}), double)
    transpose_region(Perm, Perm_T)

    --compute final coefficient matrix    
    dgemm(rank, n, n, P_final, P, Perm_T)

  else --edge case when R is full rank, P(:,I) is the identity matrix

    --transpose the Perm matrix
    var Perm_T = region(ispace(int2d, {x = n, y = n}), double)
    transpose_region(Perm, Perm_T)

    --compute final coefficient matrix
    dgemm(rank, n, n, P_final, P, Perm_T)

  end

end

--task to compute the final P matrix from P_i and P'
task id_la.compute_final_P(P_final : region(ispace(int2d), double), P_i : region(ispace(int2d), double), 
                           P_prime : region(ispace(int2d), double), N : int, target_rank : int, 
                           k : region(ispace(int2d), int))
where reads(P_i, P_prime, P_final, k), writes(P_final)
do
 
  var bounds = P_final.bounds
  var m = bounds.hi.x - bounds.lo.x + 1
  var n = bounds.hi.y - bounds.lo.y + 1
  var k_bounds = k.bounds.lo
  var final_rank = k[k_bounds]

  var P_i_bounds = P_i.bounds
  var n_P_i = P_i_bounds.hi.y - P_i_bounds.lo.y + 1
  
  --Partition P_i into individual submatrices (the diagonal submatrices of the entire matrix)
  var P_i_colors = ispace(int2d, {x = 1, y = N}) 
  var P_i_part = partition(equal, P_i, P_i_colors)

  --Partition P_prime into blocks to perform multiplication
  var P_prime_actual = region(ispace(int2d, {x = final_rank, y = N*target_rank}), double)
  for i in P_prime_actual do
    P_prime_actual[i] = P_prime[{x = i.x + P_prime.bounds.lo.x, y = i.y + P_prime.bounds.lo.y}]
  end

  var P_prime_actual_part = partition(equal, P_prime_actual, P_i_colors)

  --partition P_final for regions that are passed into DGEMM
  fill(P_final, 0.0)
  var P_final_part = partition(equal, P_final, P_i_colors) 

  --perform block matrix multiplication
  for i in P_i_colors do
    dgemm(final_rank, n_P_i/N, target_rank, P_final_part[i], P_prime_actual_part[i], P_i_part[i])
  end


end

--task to find the Frobenius norm of the difference between exact and approximate solutions
task id_la.error_frob_norm(domain_restart : region(ispace(int2d), double), domain_no_restart : region(ispace(int2d), double))
where reads(domain_restart, domain_no_restart)
do

  var rs_bounds = domain_restart.bounds
  var x_dim = rs_bounds.hi.x - rs_bounds.lo.x + 1
  var y_dim = rs_bounds.hi.y - rs_bounds.lo.y + 1

  var no_rs_bounds = domain_no_restart.bounds

  --find the difference between the two solutions
  var domain_diff = region(ispace(int2d, {x = x_dim, y = y_dim}), double)

  for i in domain_diff do
    domain_diff[i] = domain_restart[rs_bounds.lo + i] - domain_no_restart[no_rs_bounds.lo + i] 
  end

  --compute the Frobenius norm of the difference matrix
  var norm_frob : double = 0.0
 
  for i in domain_diff do
    norm_frob = norm_frob + domain_diff[i]*domain_diff[i]
  end

  --compute the Frobenius norm of the "exact" solution (used for normalization)
  var norm_frob_A : double = 0.0
 
  for i in domain_no_restart do
    norm_frob_A = norm_frob_A + domain_no_restart[i]*domain_no_restart[i]
  end

  return sqrt(norm_frob)/sqrt(norm_frob_A)

end

--task to find the maximum norm of the difference between exact and approximate solutions
task id_la.error_max_norm(domain_restart : region(ispace(int2d), double), domain_no_restart : region(ispace(int2d), double))
where reads(domain_restart, domain_no_restart)
do

  var rs_bounds = domain_restart.bounds
  var x_dim = rs_bounds.hi.x - rs_bounds.lo.x + 1
  var y_dim = rs_bounds.hi.y - rs_bounds.lo.y + 1

  var no_rs_bounds = domain_no_restart.bounds

  --find the difference between the two solutions
  var domain_diff = region(ispace(int2d, {x = x_dim, y = y_dim}), double)

  for i in domain_diff do
    domain_diff[i] = domain_restart[rs_bounds.lo + i] - domain_no_restart[no_rs_bounds.lo + i] 
  end 

  --compute the maximum norm of the difference matrix
  var norm_max : double = 0.0
 
  for i in domain_diff do
    if abs(domain_diff[i]) > norm_max then
      norm_max = abs(domain_diff[i])
    end
  end

  --compute the maximum norm of the "exact" solution (used for normalization)
  var norm_max_A : double = 0.0
 
  for i in domain_no_restart do
    if abs(domain_no_restart[i]) > norm_max_A then
      norm_max_A = abs(domain_no_restart[i])
    end
  end

  return norm_max/norm_max_A

end

return id_la
