describe "cd" do
  it "finds the correct sign-vector Z" do
    # Example from paper, page 106
    x = NMatrix[[2, -2], [0, 3], [-4, 2], dtype: :float64]
    # see printout with ssv(x, debug: true)
    assert ssv(x) == NMatrix[[-1], [1], [1], dtype: :float64]
  end

  it "correctly computes CD for the given example" do
    # Example 2 page 103
    x = NMatrix[[2, -2], [0, 3], [-4, 2], dtype: :float64]
    l_hat = NMatrix[[-2.820, 0.217], [2.278, 1.952], [4.122, -1.735],
      dtype: :float64]
    r_hat = NMatrix[[-0.651, 0.759], [0.759, 0.651]]
    assert x.approximates? l_hat.dot(r_hat.transpose), 1e-2

    l, r = cd x
    assert l.approximates? l_hat, 1e-3
    assert r.approximates? r_hat, 1e-3
    assert x.approximates? l.dot(r.transpose), 1e-15
  end

  it "correctly computes CD for a random matrix" do
    x = NMatrix.random [5,3]
    l, r = cd x
    assert x.approximates? l.dot(r.transpose), 1e-15
  end
end
