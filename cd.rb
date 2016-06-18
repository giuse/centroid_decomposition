require 'nmatrix'

# Method called at end of file for simplicity
def main
  # Example from paper, page 106
  x = NMatrix[ [2, -2], [0, 3], [-4, 2], dtype: :float64 ]
  ssv x
end

def ssv x
  z = nil

  n = 1
  pos = nil

  while true do

    # Initialize or change sign of Z
    if pos.nil?
      z = NMatrix.ones([x.rows, 1], dtype: :float64) # NMatrix[ [1], [1], [1] ]
    else
      z[pos] *= -1
    end

    # Compute S vector
    s = x.each_row.each_with_index.collect do |r, i|
      r.transpose * z[i]
    end.reduce(:+)

    # Compute V vector
    v = x.each_row.each_with_index.collect do |r, i|
      ( ((r.dot(s) * z[i]) - r.dot(r.transpose)) * z[i] ).to_a
    end.to_nm

    # Find position where V and Z have different signs and vi is largest in absolute value
    pos = nil
    val = 0
    v.zip(z).each_with_index do |(vi, zi), i|
      if vi*zi < 0 && vi.abs > val.abs
        val = vi
        pos = i
      end
    end

    # repeat until no position to change
    break if pos.nil?

    puts "n: #{n} -- s: #{s} -- v: #{v} -- z: #{z}"
    n += 1

  end

  puts "n: #{n} -- s: #{s} -- v: #{v} -- z: #{z}"

end

# Execute
main if __FILE__ == $0
