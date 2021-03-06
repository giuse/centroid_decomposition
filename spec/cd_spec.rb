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
    l, r = cd x #, debug: true
    assert x.approximates? l.dot(r.transpose), 1e-15
  end

  describe "reconstruction" do

    describe "linear regression for a linearly correlated column" do
      x = NMatrix[[3,6,9,12], dtype: :float64].transpose
      describe "when the missing values are in center column" do
        it do
          broken = NMatrix[[3,Float::NAN,Float::NAN,12], dtype: :float64].transpose
          initialize_nans broken
          assert x.approximates? broken, 1e-15
        end
      end
      describe "when the missing values are at one end of the column" do
        it "(beginning)" do
          broken = NMatrix[[Float::NAN,Float::NAN,9,12], dtype: :float64].transpose
          initialize_nans broken
          assert x.approximates? broken, 1e-15
        end
        it "(end)" do
          broken = NMatrix[[3,6,Float::NAN,Float::NAN], dtype: :float64].transpose
          initialize_nans broken
          assert x.approximates? broken, 1e-15
        end
      end
      describe "when it has only one value" do
        it do
          fixed_val = NMatrix[[3,3,3], dtype: :float64].transpose
          broken1 = NMatrix[[3,Float::NAN,Float::NAN], dtype: :float64].transpose
          broken2 = NMatrix[[Float::NAN,3,Float::NAN], dtype: :float64].transpose
          broken3 = NMatrix[[Float::NAN,Float::NAN,3], dtype: :float64].transpose
          initialize_nans broken1
          initialize_nans broken2
          initialize_nans broken3
          assert fixed_val.approximates? broken1, 1e-15
          assert fixed_val.approximates? broken2, 1e-15
          assert fixed_val.approximates? broken3, 1e-15
        end
      end
      describe "when the column is empty" do
        it "raises an exception" do
          StandardError.assert.raised? do
            broken = NMatrix[[Float::NAN,Float::NAN,Float::NAN], dtype: :float64]
            initialize_nans broken
          end
        end
      end
    end

    it "reconstructs missing values from a highly correlated matrix" do
      x = NMatrix[[3,6,9],[5,8,11],[7,10,13],[9,12,15], dtype: :float64]
      nmissing = 2
      broken = x.dup
      nmissing.times do
        r = (rand*broken.rows).to_i
        c = (rand*broken.cols).to_i
        broken[r,c] = Float::NAN # missing values
      end

      assert x.approximates? reconstruct(broken), 1e-10
    end

    # TODO: need a better way to test this!
    it "reconstructs missing values from a matrix of non-linearly correlated columns" do
      c1 = [2,4,2,16]
      c2 = c1.collect { |v| v+2 }
      c3 = c1.collect { |v| v-2 }
      x = NMatrix[c1, c2, c3, dtype: :float64].transpose
      nmissing = 2
      broken = x.dup
      broken[0,1] = Float::NAN
      broken[2,1] = Float::NAN

      # this can't be solved with linear regression
      lin = broken.dup
      missing = initialize_nans lin
      refute x.approximates? lin, 1e-3
      # but gives a rough approximation, even with little data,
      # with a lower order of magnitude of precision w.r.t. the data,
      # by using centroid decomposition reconstruction, using the
      # correlation between lines
      assert x.approximates? reconstruct(broken.dup), 1e-3
    end
    xit "reconstructs missing values from large real-world data" do
      require 'json'
      x = JSON.load(File.read 'data.json').sort_by(&:first).to_nm
      percent_missing = 0.1 #%
      nmissing = (x.cols * x.rows / 100.0 * percent_missing).round
      broken = x.dup
      nmissing.times do
        r = (rand*broken.rows).to_i
        c = (rand*broken.cols).to_i
        broken[r,c] = Float::NAN # missing values
      end

      # this should not be solvable with linear regression
      lin = broken.dup
      missing = initialize_nans lin
      refute x.approximates? lin, 1e-3
      # but gives centroid decomposition reconstruct should work
      assert x.approximates? reconstruct(broken.dup), 1e-3
    end
  end
end
