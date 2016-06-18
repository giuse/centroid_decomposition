require 'nmatrix'

# SSV algorithm straight from the pseudocode on the paper, page 105
# TODO: refactor & optimize
def ssv x, debug: false
  pos = nil  # pos = 0 in the paper
  z = NMatrix.ones([x.rows, 1], dtype: :float64)
  puts "  - x: #{x}" if debug

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
