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

def reconstruct x
  # Questions:
  # - interpolation only in column or more complex? -> linear on column
  # - Frobenius norm: can I check only between changed values? -> not yet, we'll think about that

  # Find missing values
  # - iterate rows and columns of the matrix
  # - save indices if value.is_nan? (nan_lst)

  # Initialize missing values using interpolation

  # Loop until break (the update is less than threshold)

  # Centroid decomposition
  # - just call cd on matrix

  # Dimensionality reduction
  # - set n last columns of L to zeros (n as parameter, default one, possibly two)
  # - compute L.dot(R), then get the new approximated values
  # - (optimization: only compute missing values)

  # Update missing values in x from what has been reconstructed

  # Compute Frobenius norm of previous and updated matrices
  # - euclidean distance: square all elements, sum, square root
  # - (optimization: there are few elements of difference, everything else we know is the same, we could compute the Frobenius between this subpart)

  # Stop (break) when the two matrices are less than 1e-5 away (size of update, param)


end
