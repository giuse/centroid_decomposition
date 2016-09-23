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

def frobenius x
  Math.sqrt (x**2).reduce :+
end

# Linearly interpolate two points on a new x
def linear_interpolation x1, y1, x2, y2, new_x=nil
  new_x ||= (x1+x2)/2
  m = (y2-y1)/(x2-x1)
  q = y1 - m*x1
  m*new_x + q # returns new_y
end

def reconstruct x, minimum_update_threshold=1e-5
  # Questions:
  # - interpolation only in column or more complex? -> linear on column
  # - Frobenius norm: can I check only between changed values? -> not yet, we'll think about that

  # Find missing values
  # - iterate rows and columns of the matrix
  # - save indices if value.is_nan? (nan_lst)
  missing = x.each_with_indices.collect { |v,r,c| [r,c] if v.nan? }.compact

  # Initialize missing values using interpolation
  missing.each do |r,c|
    x[r,c] = case
    when r = 0

    end
  end

  # Loop until break (the update is less than threshold)
  loop do
    # Compute Frobenius norm of previous of x
      # TODO: there are few elements of difference between the two matrices.
      # How can we compute the Frobenius only on the subparts?
    frob_old = frobenius x

    # Centroid decomposition
    l, r = cd x

    # Dimensionality reduction
    # - set n last columns of L to zeros (n as parameter, default one, possibly two)



    # - compute L.dot(R), then get the new approximated values
    reconstr = l.dot(r) # (optimization: only compute missing values)

    # Update missing values in x from what has been reconstructed
    missing.each { |r, c| x[r,c] = reconstr[r,c] }

    # Compute Frobenius norm of updated matrix
    frob_new = frobenius x

    # Stop when the frobenius of the two matrices is less than minimum_update_threshold
    break if (frob_old - frob_new).abs < minimum_update_threshold
  end

  return x
end
