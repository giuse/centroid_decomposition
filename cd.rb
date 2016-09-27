require 'nmatrix'

# SSV algorithm straight from the pseudocode on the paper, page 102
# TODO: refactor & optimize
def cd x, debug: false
  r = []
  l = []
  x.cols.times do
    z = ssv x, debug: debug
    ci = x.transpose.dot z
    ri = ci / ci.norm2
    r << ri.to_flat_a
    li = x.dot ri
    l << li.to_flat_a
    x -= li.dot ri.transpose
  end
  [l.transpose.to_nm, r.transpose.to_nm]
end

# SSV algorithm straight from the pseudocode on the paper, page 105
# TODO: refactor & optimize
def ssv x, debug: false
  pos = nil  # pos = 0 in the paper
  z = NMatrix.ones([x.rows, 1], dtype: :float64)
  puts "\n  - x: #{x}" if debug

  while true do  # do-while with internal break (/return)

    # Compute S vector
    s = x.each_row.each_with_index.collect do |r, i|
      r.transpose * z[i]
    end.reduce(:+)

    # Compute V vector
    v = x.each_row.each_with_index.collect do |r, i|
      ( ((r.dot(s) * z[i]) - r.dot(r.transpose)) * z[i] ).to_a
    end.to_nm

    # Find position where V and Z have different signs and
    # vi is largest in absolute value
    pos = nil
    val = 0
    v.zip(z).each_with_index do |(vi, zi), i|
      if vi*zi < 0 && vi.abs > val.abs  # there was an errata here in the paper
        val = vi
        pos = i
      end
    end

    puts "s: #{s} -- v: #{v} -- pos: '#{pos}' -- z: #{z}" if debug

    # repeat until no position to change
    return z if pos.nil?
    # else change sign at found position
    z[pos] *= -1
  end
end

# Frobenius norm for a matrix (generalization of L2 norm)
# @param [NMatrix]
def frobenius x
  # TODO: verify if checking only the changed values is sufficient
  Math.sqrt (x**2).reduce :+
end

# Linear interpolation function between two points
def linear_interpolation_function x1, y1, x2, y2
  m = (y2.to_f-y1)/(x2-x1)
  q = y1 - m*x1
  return -> (x) { m*x + q }
end

# Initialize contiguous missing values in a matrix column through linear interpolation
# Base case: `[...(initialized), start_value, first_nan, ..., last_nan, end_value, ...]`
# We will compute the linear regression function based on `start_value` and `end_value`.
# Then we will initialize the values in the range `first_nan..last_nan` based on that function.
def initialize_nans_block mat, col, first_nan, end_value
  # Handle border cases (->)
  # (-> empty column is handled in the caller)

  # -> column has single value, at beginning or end
  # I need to put on top a check if the single value in column is begin or end,
  if first_nan == 1 && end_value.nil?  # single value in beginning
    mat.rows.times { |row| mat[row, col] = mat[0,col] }
    return
  end
  if first_nan.zero? && end_value == mat.rows-1 # single value in end
    mat.rows.times { |row| mat[row, col] = mat[-1,col] }
    return
  end
  # I handle below the case of single-value column with value in other positions

  # -> column starts with block of nans
  if first_nan.zero?
    # as a second line-point, we need to find the next non-nan
    start_value = nil
    ((end_value+1)...(mat.rows)).each do |row|
      unless mat[row,col].nan?
        start_value = row
        break
      end
    end
    # -> only one value in the column (sub-case)
    # if we start with nan_block, find a value, then no other value found
    if start_value.nil?
      mat.rows.times { |row| mat[row, col] = mat[end_value, col] }
      return
    end
  end

  # -> column ends with block of nans
  if end_value.nil?
    # previous values are known to be already initialized
    start_value = first_nan - 1
    end_value = start_value - 1
    last_nan = mat.rows - 1
  end

  # -> base case
  # (rest of the code is shared with some of the special cases)
  start_value ||= first_nan - 1
  last_nan ||= end_value - 1
  
  # compute interpolation function based on nearest two values
  interp_fn = linear_interpolation_function(
    start_value, mat[start_value,col],
    end_value, mat[end_value,col])

  # fill the block calling the function on each NaN
  (first_nan..last_nan).each do |nan_row|
    mat[nan_row,col] = interp_fn.call(nan_row)
  end
end

# Finds and initializes missing values in a matrix
# @return [Array<Array>] missing values
def initialize_nans mat
  missing = Array.new(mat.cols) { [] }
  nan_block_begin = nil

  mat.cols.times do |col|
    mat.rows.times do |row|
      if mat[row,col].nan?
        nan_block_begin ||= row
        missing[col] << row
        # check for nan blocks at end of column
        if row==(mat.rows-1)
          # make sure the column was not entirely empty
          raise "HELL! COLUMN \##{col} IS EMPTY!" if nan_block_begin.zero?
          initialize_nans_block mat, col, nan_block_begin, nil
          nan_block_begin = nil
        end
      elsif nan_block_begin 
        # we found a value AND we're just out of a nan block
        # `nan_block_begin` points at the first `nan`
        # `r` points at the first `non-nan` (end value)
        initialize_nans_block mat, col, nan_block_begin, row
        nan_block_begin = nil
      end
    end
  end

  return missing
end

# Reconstruct missing values in `x` using Centroid Decomposition
# @param x [NMatrix] the matrix to reconstruct
# @param minimum_update_threshold termination criterion: when a reconstruction
#     iteration creates a matrix which differs from the previous by less than 
#     `minimum_update_threshold` in Frobenius norm (generalization of L2)
# @return [NMatrix] the reconstructed matrix
def reconstruct x, minimum_update_threshold=1e-5

  # Find and initialize missing values
  missing = initialize_nans x
  # Frobenius norm of "previous" x (w.r.t. update)
  frob_old = frobenius x

  # Loop until `break` (the update is less than threshold)
  loop do
    # Centroid decomposition
    l, r = cd x

    # Dimensionality reduction
    # - set n last columns of L to zeros (n as parameter, default one, possibly two)
    nzerocols = 1

    -1.downto(-nzerocols) do |col| # negative column numbering starts from last (=-1)
      l.rows.times do |row|
        l[row,col] = 0
      end
    end

    # - compute L.dot(R), then get the new approximated values
    reconstr = l.dot r.transpose # (optimization: only compute missing values?)

    # Update missing values in x from what has been reconstructed
    missing.each_with_index do |rows, c|
      rows.each do |r|
        x[r,c] = reconstr[r,c]
      end
    end

    # Compute Frobenius norm of updated matrix
    frob_new = frobenius x

    # Stop when the frobenius of the two matrices is less than minimum_update_threshold
    break if (frob_old - frob_new).abs < minimum_update_threshold

    # Update Frobenius norm of "previous" x
    frob_old = frob_new
  end

  # Finally, return reconstructed matrix
  return x
end
