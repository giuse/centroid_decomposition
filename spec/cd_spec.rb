describe "cd" do
  it "finds the correct sign-vector Z" do
    # Example from paper, page 106
    x = NMatrix[[2, -2], [0, 3], [-4, 2], dtype: :float64]
    assert ssv(x, debug: true) == NMatrix[[-1], [1], [1], dtype: :float64]
  end
end
