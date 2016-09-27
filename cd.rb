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
# Having inizialized the matrix up to this block, we expect the following sequence:
# `[start_value, first_nan, ..., last_nan, end_value]`
# We will compute the linear regression function based on `start_value` and `end_value`.
# Then we will initialize the values in the range `first_nan..last_nan` based on that function.
def initialize_nans_block mat, col, first_nan, end_value

  # TODO: handle in the caller if the column is all NANs

  # TODO: handle border cases!
  # - empty column
  # handled in the caller
  # - column starts with block of nans
  # this means that first_nan is 0, easy to check
  if first_nan.zero?
    # in that case though, we need to find the next non-nan
    # to do that, we can simply start from end_value and scroll the column until we find a number
    start_value = nil
    ((end_value+1)...(mat.rows)).each do |row|
      unless mat[row,col].nan?
        start_value = row
        break
      end
    end
    # if we end the column without finding another value, there is only one value in the column
    # we simply put the column constant with the one value we found
    if start_value.nil?
      v = mat[end_value, col]
      missing = (0...(mat.rows)).to_a - [end_value]
      missing.each { |row| mat[row,col] = v }
      return missing
    end
  end

  # - column ends with block of nans
  if end_value.nil?
    # that's easy once we're inside. All values until now are known to be initialized
    # moreover, if we're here, it means a block of nans was opened then closed, so there exist values before
    # we can just take the previous one
    start_value = first_nan - 1
    end_value = start_value - 1
    last_nan = mat.rows - 1
  end

  # with extrema taken care of, the normal case is straightforward
  start_value ||= first_nan - 1
  last_nan ||= end_value - 1

  
  interp_fn = linear_interpolation_function(
    start_value, mat[start_value,col],
    end_value, mat[end_value,col])
  (first_nan..last_nan).each do |nan_row|
    mat[nan_row,col] = interp_fn.call(nan_row)
  end
end

# Finds and initializes missing values in a matrix
# @return [Array<Array>] missing values
def initialize_nans x
  missing = Array.new(x.cols) { [] }
  nan_block_begin = nil
  # TODO: it should reset the nan_block automatically at each column end
  # the problem is most likely that each_with_indices works ROW-WISE
  # transposing x is crazy, better to use our explicit indices
  x.cols.times do |c|
    x.rows.times do |r|
      if x[r,c].nan?
  # x.each_with_indices do |v,r,c|
    # if v.nan?
      nan_block_begin ||= r
      missing[c] << r
      # HACK: refactor this
      # it's the detector for nil blocks at end of column
      if r==(x.rows-1)
        # ... which is also the detector for empty columns
        raise "HELL! EMPTY COLUMN!" if nan_block_begin==0
        initialize_nans_block x, c, nan_block_begin, nil
        nan_block_begin = nil
      end
    elsif nan_block_begin 
      # we found a value AND we're just out of a nan block
      # `nan_block_begin` points at the first `nan`
      # `r` points at the first `non-nan` (end value)
      initialize_nans_block x, c, nan_block_begin, r
      nan_block_begin = nil
    end
  end end
  return missing
end

# Reconstruct missing values in `x` using Centroid Decomposition
# @param x [NMatrix] the matrix to reconstruct
# @param minimum_update_threshold termination criterion: when a reconstruction
#     iteration creates a matrix which differs from the previous by less than 
#     `minimum_update_threshold` in Frobenius norm (generalization of L2)
# @return [NMatrix] the reconstructed matrix
def reconstruct x, minimum_update_threshold=1e-5
  # NOTE: I still need a case/switch to correctly call the interpolation:
  # - given a missing value, search for precedent and subsequent non-nil values
  # - if beginning/end of list is reached, use two subsequent/precedent
  # - nil value clusters should be updated together from the same function
  # - this means I should return the interpolating function rather than a y value
  # - last, I should go at this top-down per each column, thus maintaining a precedent/subsequent, rather than collecting the missing and reconstructing them singularly as I am doing here

  # Find and initialize missing values
  # - iterate rows and columns of the matrix -- maintain last_value
  # - if `value.nan?` save indices to missing list and initialize

  missing = initialize_nans x
  frob_old = frobenius x  # Frobenius norm of "previous" x

  # Loop until `break` (the update is less than threshold)
  loop do
    # Centroid decomposition
    l, r = cd x

    # Dimensionality reduction
    # - set n last columns of L to zeros (n as parameter, default one, possibly two)

    # TODO!

    # - compute L.dot(R), then get the new approximated values
    reconstr = l.dot(r) # (optimization: only compute missing values?)

    # Update missing values in x from what has been reconstructed
    missing.each { |r, c| x[r,c] = reconstr[r,c] }

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
