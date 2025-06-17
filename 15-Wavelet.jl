### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 89a01055-d24b-488c-a3c8-9c8822567cd8
using Plots

# ╔═╡ 359483d5-0d23-4d47-ab24-cbfaa505b16e
md"""
# Wavelet Transform

The wavelet transform works for many signals as a sparsity-promoting
transform. Similarly to the Fourier transform, the wavelet transform
can be continuous or discrete. Instead of the notion of frequency in
the Fourier transform, the wavelet transform has the notion of scale,
which makes it well suited for representing non-stationary signals
with scale-dependent behavior.
"""

# ╔═╡ e451844b-0017-4c03-95e1-561795859321
md"""
## Digital wavelet transform

Similarly to the Fourier transform, there is a continuous analog of the digital wavelet transform, known as the continuous wavelet transform. In this chapter, we will limit the discussion to the digital version.
"""

# ╔═╡ a9957069-99bd-4c4a-8e31-72e3f3bbf87a
md"""
The digital wavelet transform (DWT) is usually not derived directly from the continuous wavelet transform but is constructed using operations of special digital filters. One of the simplest versions of DWT is the Haar transform.
"""

# ╔═╡ 28a36422-722f-4726-954b-a68b40b58705
md"""
## Haar transform

Let us consider a simple example to see how the Haar transform works. Suppose that the input data contains the following set of numbers (summer weather forecast with temperature in Fahrenheit):

$$\begin{array}{|cccccccc|} \hline 94 & 94 & 95 & 97 & 98 & 98 & 96 & 96 \\ \hline \end{array}\;.$$
"""

# ╔═╡ 939e0835-bafe-4ae1-be4f-519a2d030580
md"""
The first step of the Haar transform is to divide the numbers into pairs:

$$\begin{array}{|cc|cc|cc|cc|} \hline 94 & 94 & 95 & 97 & 98 & 98 & 96 & 96 \\ \hline \end{array}\;.$$

Next, let us compute averages and differences between the two numbers in each pair:

$$\begin{array}{|cc|cc|cc|cc|} \hline 94 & 0 & 96 & 2 & 98 & 0 & 96 & 0 \\ \hline \end{array}\;.$$
"""

# ╔═╡ 202fcfd9-cb0d-4153-9672-bdfc008c435e
md"""
Now we are going up the scale ladder: let us put the differences
aside, concentrate on the averages, and consider them as the new input:
Next, let us compute averages and differences between the two numbers in each pair:

$$\begin{array}{|cccc|cccc|} \hline 94 & 96 & 98 & 96 & 0 & 2 & 0 & 0 \\ \hline \end{array}\;.$$

Again, divide into pairs:

$$\begin{array}{|cc|cc|cccc|} \hline 94 & 96 & 98 & 96 & 0 & 2 & 0 & 0 \\ \hline \end{array}\;.$$
"""

# ╔═╡ 80cc744b-3209-4129-9eee-b9239b6cd601
md"""
Compute averages and differences:

$$\begin{array}{|cc|cc|cccc|} \hline  95 & 2 & 97 & -2 & 0 & 2 & 0 & 0 \\ \hline \end{array}\;.$$

Put differences aside and concentrate on the averages:

$$\begin{array}{|cc|cccccc|} \hline  95 & 97 & 2 & -2 & 0 & 2 & 0 & 0 \\ \hline \end{array}\;.$$
"""

# ╔═╡ 8b952246-80b2-44f6-b8d2-c8d6f80adfda
md"""
Now we have just one pair. The process ends with computing the last average and the last difference:

$$\begin{array}{|c|ccccccc|} \hline 96 & 2 & 2 & -2 & 0 & 2 & 0 & 0 \\ \hline \end{array}\;.$$
"""

# ╔═╡ 5a26a2da-030b-458f-b123-36fee9b13ba7
md"""
The first coefficient in the transformed data represents the global average, while the other coefficients stand for differences at different scales. Because the input data are smooth at different scales, the differences are minor in size and take a small dynamic range, which makes them easily compressible. The Haar transform of piecewise smooth data is sparse in the sense of having only a small number of significant coefficients.
"""

# ╔═╡ 7116a716-51f0-4849-810a-9b834d8e31cf
function haar!(x::Vector{T}, inv::Bool=false) where T <: Real
    n = length(x) # Assume n is a power of 2
    if inv
        j = n ÷ 2
        while j >= 1
            for i in 1:2j:n-j
                x[i] -= x[i+j]/2
                x[i+j] += x[i]
            end
            j ÷= 2
        end
    else
        j = 1
        while j <= n ÷ 2
            for i in 1:2j:n-j
                x[i+j] -= x[i]
                x[i] += x[i+j]/2
            end
            j *= 2
        end
    end
end

# ╔═╡ 291f65a8-a649-42f4-99ca-2fa3c66ecef0
begin
	x = [94, 94, 95, 97, 98, 98, 96, 96]
	haar!(x)
	x
end

# ╔═╡ 0e81369e-bd42-411e-b767-113097242f8c
begin
	y = [94, 94, 95, 97, 98, 98, 96, 96]
	haar!(y)
	haar!(y, true)
	y
end

# ╔═╡ 5ce31f40-776e-49fc-ae87-6d8b908c645a
function order_coefficients(c::Vector{T}, n::Int) where T <: Real
    "put wavelet coefficients in the right order"
    nc = length(c)
    x = Array{T}(undef, n)
    ic = 1
    x[1] = c[1]
    j = nc ÷ 2
    while j >= 1                  
        for i in 1:2j:nc-j
            ic += 1
            x[ic] = c[i+j]
        end
        if ic >= n
            return x
        end
        j ÷= 2
    end
    return x
end

# ╔═╡ 57fc5904-b45b-4390-a250-7a9cff4c8c5a
function reorder_coefficients(x::Vector{T}, nc::Int) where T <: Real
    "inverse transformation to order_coefficients"
    n = length(x)
    c = Array{T}(undef, nc)
    ic = 1
    c[1] = x[1]
    j = nc ÷ 2
    while j >= 1
        for i in 1:2j:nc-j
            if ic < n
                ic += 1
                c[i+j] = x[ic]
            else
                c[i+j] = zero(T)
            end
        end
        j ÷= 2
    end  
    return c
end

# ╔═╡ 352401a1-9016-407e-88ee-2dfe784f9057
function DWT(transform::Function, x::Vector{T}, 
             inv::Bool=false) where T <: Real
    y = similar(x)
    n = length(x)
    # find the nearest power of two
    nt = 1
    while nt < n
        nt *= 2
    end
    if inv
        t = reorder_coefficients(x, nt)
    else
        t = vcat(x, zeros(T, nt-n)) # pad with zeros
    end
    transform(t, inv)     
    if inv
        y = t[1:n] # truncate
    else  
        y = order_coefficients(t, n)
    end   
    return y
end

# ╔═╡ 2de4580d-319e-4b96-bcee-c4dc7e65aaad
DWT(haar!, [94, 94, 95, 97, 98, 98, 96, 96])

# ╔═╡ 4abb361a-d914-4fc3-a422-4a0341f1e4e0
DWT(haar!, [96, 2, 2, -2, 0, 2, 0, 0], true)

# ╔═╡ 3b655702-ffb6-4d0a-8360-02d94e1f080d
download("https://ahay.org/data/imgs/mona.img","mona.img")

# ╔═╡ d7a8135b-df9d-4882-b277-7489d1229b6a
begin
	image = Array{UInt8}(undef, 512, 513) # byte array
	read!("mona.img", image)
	image = reverse(transpose(image[:,1:512]), dims=1)
end

# ╔═╡ c495f591-f815-4189-bd0c-82fa209a5fff
heatmap(image, legend=:none, box=:none, 
        aspect_ratio=:equal, xlim=[0.5,512.5])

# ╔═╡ 230dd685-c9fe-4785-8eb6-177ad8d9a521
Haar(trace) = DWT(haar!, trace)

# ╔═╡ cfb14a1b-7274-467b-b588-493f0069c834
haar1 = mapslices(Haar, float.(image); dims=2)

# ╔═╡ f0a1d712-14bc-47fc-894c-8e9f0525c9c1
function plot_wavelet1(dwt, name)
    plt = heatmap(dwt, legend=:none, box=:none, color=:grays, 
                  aspect_ratio=:equal, xlim=[0.5,512.5],
                  title="1-D $name Transform", clim=(0, 42))
    for s in (256, 128, 64, 32, 16, 8, 4, 2)
        plot!(plt, [s, s], [0, 512], legend=:none, color=:yellow)
    end
    return plt
end

# ╔═╡ 9a69f3f1-b85e-4218-a06b-ee290d41f3a4
plot_wavelet1(haar1, "Haar")

# ╔═╡ 4b6430f7-25eb-4d8b-8b77-1e3a83f58f37
haar2 = mapslices(Haar, haar1; dims=1)

# ╔═╡ 2baf629c-df91-45cc-a994-427b8d64d8a1
function plot_wavelet2(dwt, name)
    plt = heatmap(dwt, legend=:none, box=:none, color=:grays, 
                  aspect_ratio=:equal, xlim=[0.5,512.5],
                  title="2-D $name Transform", clim=(0, 42))
    for s in (256, 128, 64, 32, 16, 8, 4, 2)
        plot!(plt, [s, s], [0, 512], legend=:none, color=:yellow)
        plot!(plt, [0, 512], [s, s], legend=:none, color=:yellow)
    end
    return plt
end

# ╔═╡ 47f54019-42b4-4598-84e0-1a56bc988de1
plot_wavelet2(haar2, "Haar")

# ╔═╡ 240a7adb-bb7e-4919-8ee5-e98e1c618482
md"""
Before we move to a more general formulation, let us make some
observations about the properties of the transform:

1. The cost of the operation satisfies the recursion

    $$C(N) = C(N/2) + \gamma\,N\;,$$

  which ends up in the cost $O(N)$, as we saw before. Thus, the digital wavelet transform is faster than the fast Fourier transform, which costs $O(N\,\log\,N)$.
"""

# ╔═╡ b859db26-6a71-4cce-9ea0-19f738680c71
md"""
2. If the digital Fourier transform decomposes signals by frequencies, the leading characteristic of the wavelet transform is scale. 
"""

# ╔═╡ ba630263-19b1-47f1-b6ca-ee5843c49df8
md"""

3. The operation is easily invertible. From the average and the difference between the odd and even values

$$\begin{array}{rcl}
    a & = & (o+e)/2\;, \\
    d & = & o-e\;,\end{array}$$

it is easy to reconstruct both values

$$\begin{array}{rcl}
    o & = & a+d/2\;. \\
    e & = & a - d/2\;.\end{array}$$
"""

# ╔═╡ 748ff555-8c4e-415e-bf6d-33be6c49591b
md"""
If the forward operation moves from small scale to large scale, the inverse transform moves in the opposite direction, reconstructing odd and even values at different scales going from large to small.      
"""

# ╔═╡ eb22ac80-cee8-47b3-9057-f7c1bb4d7a1f
#inverse transform
Haari(trace) = DWT(haar!, trace, true)

# ╔═╡ af89270e-bde4-4daa-be11-96b72b61f484
inverse = mapslices(Haari, mapslices(Haari, haar2; dims=1); dims=2);

# ╔═╡ 5772dd5b-b4b1-4bbc-86ce-97dededbd613
heatmap(inverse, legend=:none, box=:none, aspect_ratio=:equal, 
        xlim=[0.5,512.5], title="Inverse Transform")

# ╔═╡ 2448b7b0-d7e4-4950-9e53-44a89e7d41db
md"""
## Filtering formulation

If a digital signal is represented by its Z-transform

$$X(Z)=x_0+x_1\,Z+x_2\,Z^2+\cdots \;,$$

the operation of upsampling (inserting zeros between each two samples)
can be represented by

$$X_{\uparrow}(Z)=X\left(Z^2\right) = x_0+x_1\,Z^2+x_2\,Z^4+\cdots$$
"""

# ╔═╡ c9037f5a-6847-4d1e-8d9d-3b8b94aa522d
md"""
The inverse operation is subsampling (selecting every other sample.) Noting that

$$X(-Z)=x_0-x_1\,Z+x_2\,Z^2+\cdots \;,$$

we can find that

$$\begin{array}{rcl}X_{\downarrow}(Z) & = & \displaystyle \frac{1}{2}\left[X\left(\sqrt{Z}\right)+X\left(-\sqrt{Z}\right)\right] \\
   & = & \displaystyle x_0+x_2\,Z+x_4\,Z^2+\cdots\end{array}$$
"""

# ╔═╡ dcab260a-cd2d-4603-9ceb-787a50aff76d
md"""
The operation of converting the signal into odd and even components and then performing filtering operations (averages and differences in the case of the Haar transform) can be described as a chain of filtering by two filters: low-pass $L(Z)$ and high-pass $H(Z)$, followed by subsampling. 
"""

# ╔═╡ 340467eb-0b82-43e9-b0e4-429cdc6e6064
md"""
Thus, the output from one step of DWT is

$$\begin{array}{rcl} A(Z) & = & \displaystyle
   \frac{1}{2}\left[L\left(\sqrt{Z}\right)\,X\left(\sqrt{Z}\right)\right. \\ & & \displaystyle \left. + L\left(-\sqrt{Z}\right) X\left(-\sqrt{Z}\right)\right]\;,
   \\
   D(Z) & = & \displaystyle
   \frac{1}{2}\left[H\left(\sqrt{Z}\right)\,X\left(\sqrt{Z}\right)\right. \\ & & \displaystyle \left. +H\left(-\sqrt{Z}\right)
              X\left(-\sqrt{Z}\right)\right]\;.\end{array}$$

At the next step, $X(Z)$ is replaced by (twice smaller) $A(Z)$, and the operation is repeated.
"""

# ╔═╡ 96f1c5f7-f443-4095-a222-e80451588ab2
md"""
The inverse transform will involve upsampling followed by additional
filtering by filters $\hat{L}(Z)$ and $\hat{H}(Z)$ and summing the
results. The output takes the form

$$\begin{array}{rcl}
  \hat{X}(Z) & = & \hat{H}(Z)\,A(Z^2) + \hat{L}(Z)\,D(Z^2) \\ & = & \displaystyle
  X(Z)\,\frac{\hat{H}(Z)\,L(Z)+\hat{L}(Z)\,H(Z)}{2} \\ & & \displaystyle
  + X(-Z)\,\frac{\hat{H}(Z)\,L(-Z)+\hat{L}(Z)\,H(-Z)}{2}\;.\end{array}$$
"""

# ╔═╡ 480411a7-aab9-4837-bd96-36cfd0d1992a
md"""
To guarantee that the signal is reconstructed exactly $\hat{X}(Z) =
X(Z)$, we need to find a pair of *analysis* filters $L(Z)$ and $H(Z)$ and a pair
of *synthesis* filters $\hat{L}(Z)$ and $\hat{H}(Z)$that satisfy two conditions:

$$\begin{array}{rcl}
  \hat{H}(Z)\,L(Z)+\hat{L}(Z)\,H(Z) & = & 2\;, \\
  \hat{H}(Z)\,L(-Z)+\hat{L}(Z)\,H(-Z) & = & 0\;.\end{array}$$
"""

# ╔═╡ 9556acd8-154e-4dc5-8d6e-418992d5ee0b
md"""
We can satisfy the second condition by choosing

$$\begin{array}{rcl}
  \hat{L}(Z) & = & -L(-Z)\,Z^{-k}\;, \\
   \hat{H}(Z) & = & H(-Z)\,Z^{-k}\;,\end{array}$$       

where $k$ is some integer power. This choice leaves one condition on the
connection between $L(Z)$ and $H(Z)$, which takes the form

$$L(Z)\, H(-Z) - L(-Z)\,H(Z) = 2\,Z^{k}\;.$$
"""

# ╔═╡ 24f792a7-14ad-4ca7-8691-7d601eb380e3
md"""
In the case of the Haar transform, 

$$\begin{array}{rcl}  L(Z) & = & \displaystyle \frac{1+Z}{2} \\ H(Z) & = & 1-Z\;,\end{array}$$ 

and it is easy to check that the condition is satisfied with $k=1$.
"""

# ╔═╡ aa5770df-7711-447a-85f7-832133c1ab73
md"""
Another example is

$$\begin{array}{rcl}  L(Z) & = & \displaystyle \frac{(1+Z)^2\,(-1+4Z-Z^2)}{4}\;, \\
  H(Z) & = & \displaystyle \frac{(1-Z)^2}{4}\;,\end{array}$$

which satisfies the invertibility condition with $k=3$.
"""

# ╔═╡ 5f3dae3a-202d-4378-b641-6f2ae5266d7d
md"""
An additional symmetry can be imposed by taking $H(Z)=L(-Z)$ so that the invertibility condition reduces to

$$L^2(Z) - L^2(-Z) = 2\,Z^{k}\;.$$

Filters that satisfy this condition are known as *quadrature-mirror filters*.The Haar filters can satisfy it if defined symmetrically with$L(Z) = (1+Z)/\sqrt{2}$ and $H(Z) = (1-Z)/\sqrt{2}$. 
"""

# ╔═╡ 69802871-9d63-4193-a392-d26a89c8724f
md"""
The cost of DWT grows proportionally to the length of the filter. On the other hand, longer filters allow us to achieve better signal compression while retaining the main properties of the transform: $O(N)$ cost, scale separation, and exact invertibility.
"""

# ╔═╡ bbeb09ca-af66-4897-a780-0adc0a99a680
function wavelet!(x::Vector{T}, inv::Bool=false) where T <: Real
    n = length(x) # Assume n is a power of 2
    if inv
        j = n ÷ 2
        while j >= 1
            for i in 1+2j:2j:n-j; x[i] -= (x[i + j] + x[i - j])/4; end
            x[1] -= x[j]/2
            for i in 1:2j:n-2j; x[i + j] += (x[i] + x[i + 2j])/2; end
            if n > 2j; x[n - j] += x[n - 2j]; end
            j ÷= 2
        end
    else
        j = 1
        while j <= n ÷ 2 
            if n > 2j; x[n - j] -= x[n - 2j]; end
            for i in 1:2j:n-2j; x[i+j] -= (x[i] + x[i + 2j])/2; end
            x[1] += x[j]/2
            for i in 1+2j:2j:n-j; x[i] += (x[i + j] + x[i - j])/4; end
            j *= 2
        end
    end
end

# ╔═╡ 1b7fd166-bd10-4af2-a9a4-7ddb0f86d2a1
Wavelet(trace) = DWT(wavelet!, trace)

# ╔═╡ ea3ed115-dc10-4754-971e-0b80f8608889
dwt1 = mapslices(Wavelet, float.(image); dims=2);

# ╔═╡ d5606f38-9231-463b-93fa-b74a7651fbcf
dwt2 = mapslices(Wavelet, dwt1; dims=1);

# ╔═╡ a76601ad-6e79-4b07-9e99-5285a5e090db
plot_wavelet1(dwt1, "Wavelet")

# ╔═╡ 6dcd45b7-fdd4-4e7e-8ec6-e1efae25c78d
plot_wavelet2(dwt2, "Wavelet")

# ╔═╡ b201e2fe-e14a-4d76-8266-181a40b3c89a
coef(transform) = sort(abs.(vec(transform)), rev=true)

# ╔═╡ 9a9d6b70-8914-4aa1-bfde-57620f25ce23
begin
	coef1=coef(haar2)
	coef2=coef(dwt2)
end

# ╔═╡ 078571af-5e5b-4ba2-86b5-50949f16b737
plot([coef1[1:5000] coef2[1:5000]], labels=["Haar" "Wavelet"], 
     linewidth=2, yscale=:log2, title="Decay of Coefficients")

# ╔═╡ eefd1fc3-c0e3-4cfb-b42b-fb2ead981b33
md"""
## Lifting scheme

The high-pass filter in the Haar transform $H(Z)=1-Z$ can be interpreted as *prediction and subtraction*. Each sample is predicted by its neighboring sample, and the filter outputs the prediction residual. Similarly, we can rewrite the high-pass filter $H(Z)$ in the other version of DWT as a scaled and shifted version of $1-(Z^{-1}+Z)/2$, where we predict each sample by a mean average of its left and right neighbors. 
"""

# ╔═╡ 10a9dbab-eb69-4489-9f50-18a915468a78
md"""
This observation suggests an idea or writing DWT in the following form:

$$\begin{array}{rcl}\mathbf{r}  & = & \mathbf{o} - \mathbf{P[e]}\; \\
    \mathbf{a}  & = & \mathbf{e} + \mathbf{U[r]}\;,\end{array}$$

where $\mathbf{o}$ and $\mathbf{e}$ represent the odd an even components of the signal, $\mathbf{P}$ and $\mathbf{U}$ are *prediction* and *update* convolutional operators, and  $\mathbf{d}$ and $\mathbf{a}$ are output differences and averages. At the next scale, the average $\mathbf{a}$ becomes the new signal and is split into differences and averages again.
"""

# ╔═╡ ce7de978-7e46-4851-a94f-8f0be028a438
md"""
This rearrangement, known as the *lifting scheme*, allows the computation to proceed *in place* without creating new copies of the data. The data vector is replaced by its transform without duplication. It is also easy to write the inverse operation by simply reversing the steps as follows:

$$\begin{array}{rcl}\mathbf{e} & = & \mathbf{a} - \mathbf{U[r]}\;, \\
    \mathbf{o} & = & \mathbf{r}  + \mathbf{P[e]}\;,\end{array}$$

To define the transform, we need to choose an appropriate $\mathbf{P}$ and design $\mathbf{U}$ to satisfy the previously defined condition.
"""

# ╔═╡ e675f460-d681-40b4-8479-2c722d7a9b5f
md"""
As an example, the filters we used in the example above use the following forms of prediction and update:

$$\begin{array}{rcl}  \mathbf{P[e]}_n & = & \left(\mathbf{e}_{n-1} +
    \mathbf{e}_{n}\right)/2\;, \\ 
  \mathbf{U[r]}_n & = & \left(\mathbf{r}_{n-1} + \mathbf{r}_{n}\right)/4\;.,\end{array}$$

Longer DWT filters may require multiple lifting steps.
"""

# ╔═╡ 5e4a93ff-fe0a-4dda-a979-ee6add0ba3d7
md"""
## From wavelets to seislets

The idea of the seislet transform is to use the template of DWT but to reimagine the prediction and update steps in the lifting scheme so
that they would correspond to predicting the most essential features in the signal.
"""

# ╔═╡ 6a618dfa-5123-47bf-9ac2-5d7fc5c32fe1
md"""
Note, for example, that the filter $H(Z)=1-Z$ in the Haar transform corresponds to the exact prediction of zero-frequency signals. To modify the prediction for a non-zero frequency $\omega_0$, it is sufficient to modify the prediction by modulation $H(Z)=1-Z/Z_0$, where $Z_0=e^{-i\,\omega_0\,\Delta t}$. Changing $Z$ to $Z/Z_0$ in the definition of DWT defines the one-dimensional seislet transform.
"""

# ╔═╡ 9df0a409-1c37-4ea9-b935-f0b19ef842d9
md"""
What if we want to use multiple frequencies instead of just one dominant frequency? Using a seislet transform for each frequency generates an overcomplete representation, known as a *frame*. Overcompleteness allows us to search for the sparsest possible
representation of the data in the transform domain.
"""

# ╔═╡ a9d23ab0-216b-4d5e-82f0-7933d5dcb552
md"""
## References

The idea of the Haar transform was proposed in 1909 by Alfréd Haar, a Hungarian mathematician. It took a long time before more general concepts in wavelet analysis could be developed and applied to digital data.
"""

# ╔═╡ d0cef261-f6f0-4692-b217-e7e7ce29644f
md"""
Morlet (1981) first proposed the continuous wavelet transform, which was motivated, together with the term "wavelet," by seismic data analysis. Some of the world's brightest applied mathematicians contributed to developing a multiresolution theory that bridges the gap between continuous wavelets and filtering techniques from digital signal processing (Daubechies, 1992). For a comprehensive treatment of the wavelet theory, see (Mallat, 2008). Jensen and la Cour-Harbo (2001) provide a concise introduction to DWT. 
"""

# ╔═╡ e8955866-0d1c-4423-89c6-cb2c521b3367
md"""
* Daubechies, I., 1992, Ten lectures on wavelets: SIAM.
* Jensen, A., and A. la Cour-Harbo, 2001, Ripples in mathematics: the discrete wavelet transform: Springer Science & Business Media.
* Mallat, S., 2008, A wavelet tour of signal processing: The sparse way, 3rd ed.: Academic Press, Inc.
* Morlet, J., 1981, Sampling theory and wave propagation: 51st Ann. Internat. Mtg, Soc. of Expl. Geophys., Session:S15.1.
"""

# ╔═╡ 0c551855-1ca7-4713-a916-cf97f590157d
md"""

The lifting scheme was proposed by Sweldens (1995) and further developed by Daubechies & Sweldens (1998).

* Daubechies, I., and W. Sweldens, 1998, Factoring wavelet transforms into lifting steps: Journal of Fourier analysis and applications, 4, 247–269.
* Sweldens, W., 1995, The lifting scheme: A new philosophy in biorthogonal wavelet con- structions: Wavelet Applications in Signal and Image Processing III, Proc. SPIE 2569, 68–79.
"""

# ╔═╡ db5d1895-19e5-4457-b999-37a10c4fb708
md"""
The field of image analysis brought several wavelet-like representations that explore the directional characteristics of images (Welland, 2003). Among those transforms are contourlets (Do and Vetterli, 2005), shearlets (Guo & Labate, 2007), and some others. These transforms attempt to represent elongated image features, such as curved edges.  
"""

# ╔═╡ 75e78f81-8fe4-4c85-b335-029a28176f79
md"""
* Do, M. N., and M. Vetterli, 2005, The contourlet transform: an efficient directional multiresolution image representation: IEEE Transactions on Image Processing, 14, 2091–2106. 
* Guo, K., and D. Labate, 2007, Optimally sparse multidimensional representation using shearlets: SIAM Journal on Mathematical Analysis, 39, 298–318.
* Welland, G., ed., 2003, Beyond wavelets: Academic Press.
"""

# ╔═╡ 568df3d4-16a2-42d7-bae3-c5ef69a0d967
md"""
Curvelets (Starck et al., 2002; Candès et al., 2006) have become especially popular for analyzing seismic images (Herrmann & Hennenfent, 2008).
"""

# ╔═╡ d40fb873-15c1-4ead-9bf7-9a9df9cb6190
md"""
* Candès, E., Demanet, L., Donoho, D. and Ying, L., 2006. Fast discrete curvelet transforms. multiscale modeling & simulation, 5, 861-899.
* Herrmann, F. J., and G. Hennenfent, 2008, Non-parametric seismic data recovery with curvelet frames: Geophysical Journal International, 173, 233–248.
* Starck, J. L., E. J. Candès, and D. L. Donoho, 2002, The curvelet transform for image denoising: IEEE Transactions on Image Processing, 11, 670–684.
"""

# ╔═╡ 4a51e4b2-e388-4c8c-9e2a-b3d46fc6f651
md"""
The seislet transform was proposed by Fomel (2006) and extended by (Fomel & Liu, 2010).

* Fomel, S., 2006, Towards the seislet transform, in 76th Ann. Internat. Mtg: Soc. of Expl. Geophys., 2847–2850.
* Fomel, S., and Y. Liu, 2010, Seislet transform and seislet frame: Geophysics, 75, V25–V38. 
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
Plots = "~1.40.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.4"
manifest_format = "2.0"
project_hash = "433c825682ac2300f10e88ef6fc0d59d12a6d791"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8873e196c2eb87962a2048b3b8e08946535864a1"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+4"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "c785dfb1b3bfddd1da557e861b919819b82bbe5b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.27.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "64e15186f0aa277e174aa81798f7eb8598e0157e"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "f36e5e8fdffcb5646ea5da81495a5a7566005127"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.3"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e51db81749b0777b2147fbe7b783ee79045b8e99"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.4+3"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "21fac3c77d7b5a9fc03b0ec503aa1a6392c34d2b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.15.0+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "786e968a8d2fb167f2e4880baba62e0e26bd8e4e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.3+1"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "52adc6828958ea8a0cf923d53aa10773dbca7d5f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.9"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4e9e2966af45b06f24fd952285841428f1d6e858"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.9+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "48b5d4c75b2c9078ead62e345966fa51a25c05ad"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.2+1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "01979f9b37367603e2848ea225918a3b3861b606"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "627fcacdb7cb51dc67f557e1598cdffe4dda386d"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.14"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "55c53be97790242c29031e5cd45e8ac296dadda3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "71b48d857e86bf7a1838c4736545699974ce79a2"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.9"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "be3dc50a92e5a386872a493a10050136d4703f9b"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "36bdbc52f13a7d1dcb0f3cd694e01677a515655b"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "854a9c268c43b77b0a27f22d7fab8d33cdb3a731"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+3"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "8be878062e0ffa2c3f67bb58a595375eda5de80b"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.11.0+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "ff3b4b9d35de638936a525ecd36e86a8bb919d11"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c6ce1e19f3aec9b59186bdf06cdf3c4fc5f5f3e6"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.50.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "61dfdba58e585066d8bce214c5a51eaa0539f269"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "84eef7acd508ee5b3e956a2ae51b05024181dee0"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.2+2"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "b404131d06f7886402758c9ce2214b636eb4d54a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "edbf5309f9ddf1cab25afc344b1e8150b7c832f9"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.2+2"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+4"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+3"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "12f1439c4f986bb868acda6ea33ebc78e19b95ad"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.7.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "dae01f8c2e069a683d3a6e17bbae5070ab94786f"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.9"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "01915bfcd62be15329c9a07235447a89d588327c"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.1"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "a2fccc6559132927d4c5dc183e3e01048c6dcbd6"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "7d1671acbe47ac88e981868a078bd6b4e27c5191"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.42+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "15e637a697345f6743674f1322beefbc5dcd5cfc"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.3+2"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "9dafcee1d24c4f024e7edc92603cedba72118283"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+3"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2b0e27d52ec9d8d483e2ca0b72b3cb1a8df5c27a"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+3"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "02054ee01980c90297412e4c809c8694d7323af3"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+3"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d7155fea91a4123ef59f42c4afb5ab3b4ca95058"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+3"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a1a7eaf6c3b5b05cb903e35e8372049b107ac729"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.5+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "b6f664b7b2f6a39689d822a6300b14df4668f0f4"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.4+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a490c6212a0e90d2d55111ac956f7c4fa9c277a6"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+1"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee57a273563e273f0f53275101cd41a8153517a"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "1a74296303b6524a0472a8cb12d3d87a78eb3612"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "dbc53e4cf7701c6c7047c51e17d6e64df55dca94"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+1"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "ab2221d309eda71020cdda67a973aa582aa85d69"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+1"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b9ead2d2bdb27330545eb14234a2e300da61232e"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+3"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6e50f145003024df4f5cb96c7fce79466741d601"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.56.3+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0ba42241cb6809f1a278d0bcb976e0483c3f1f2d"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+1"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "b70c870239dc3d7bc094eb2d6be9b73d27bef280"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.44+2"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─359483d5-0d23-4d47-ab24-cbfaa505b16e
# ╟─e451844b-0017-4c03-95e1-561795859321
# ╟─a9957069-99bd-4c4a-8e31-72e3f3bbf87a
# ╟─28a36422-722f-4726-954b-a68b40b58705
# ╟─939e0835-bafe-4ae1-be4f-519a2d030580
# ╟─202fcfd9-cb0d-4153-9672-bdfc008c435e
# ╟─80cc744b-3209-4129-9eee-b9239b6cd601
# ╟─8b952246-80b2-44f6-b8d2-c8d6f80adfda
# ╟─5a26a2da-030b-458f-b123-36fee9b13ba7
# ╠═7116a716-51f0-4849-810a-9b834d8e31cf
# ╠═291f65a8-a649-42f4-99ca-2fa3c66ecef0
# ╠═0e81369e-bd42-411e-b767-113097242f8c
# ╠═5ce31f40-776e-49fc-ae87-6d8b908c645a
# ╠═57fc5904-b45b-4390-a250-7a9cff4c8c5a
# ╠═352401a1-9016-407e-88ee-2dfe784f9057
# ╠═2de4580d-319e-4b96-bcee-c4dc7e65aaad
# ╠═4abb361a-d914-4fc3-a422-4a0341f1e4e0
# ╠═3b655702-ffb6-4d0a-8360-02d94e1f080d
# ╠═d7a8135b-df9d-4882-b277-7489d1229b6a
# ╠═89a01055-d24b-488c-a3c8-9c8822567cd8
# ╠═c495f591-f815-4189-bd0c-82fa209a5fff
# ╠═230dd685-c9fe-4785-8eb6-177ad8d9a521
# ╠═cfb14a1b-7274-467b-b588-493f0069c834
# ╠═f0a1d712-14bc-47fc-894c-8e9f0525c9c1
# ╠═9a69f3f1-b85e-4218-a06b-ee290d41f3a4
# ╠═4b6430f7-25eb-4d8b-8b77-1e3a83f58f37
# ╠═2baf629c-df91-45cc-a994-427b8d64d8a1
# ╠═47f54019-42b4-4598-84e0-1a56bc988de1
# ╟─240a7adb-bb7e-4919-8ee5-e98e1c618482
# ╟─b859db26-6a71-4cce-9ea0-19f738680c71
# ╟─ba630263-19b1-47f1-b6ca-ee5843c49df8
# ╟─748ff555-8c4e-415e-bf6d-33be6c49591b
# ╠═eb22ac80-cee8-47b3-9057-f7c1bb4d7a1f
# ╠═af89270e-bde4-4daa-be11-96b72b61f484
# ╠═5772dd5b-b4b1-4bbc-86ce-97dededbd613
# ╟─2448b7b0-d7e4-4950-9e53-44a89e7d41db
# ╟─c9037f5a-6847-4d1e-8d9d-3b8b94aa522d
# ╟─dcab260a-cd2d-4603-9ceb-787a50aff76d
# ╟─340467eb-0b82-43e9-b0e4-429cdc6e6064
# ╟─96f1c5f7-f443-4095-a222-e80451588ab2
# ╟─480411a7-aab9-4837-bd96-36cfd0d1992a
# ╟─9556acd8-154e-4dc5-8d6e-418992d5ee0b
# ╟─24f792a7-14ad-4ca7-8691-7d601eb380e3
# ╟─aa5770df-7711-447a-85f7-832133c1ab73
# ╟─5f3dae3a-202d-4378-b641-6f2ae5266d7d
# ╟─69802871-9d63-4193-a392-d26a89c8724f
# ╠═bbeb09ca-af66-4897-a780-0adc0a99a680
# ╠═1b7fd166-bd10-4af2-a9a4-7ddb0f86d2a1
# ╠═ea3ed115-dc10-4754-971e-0b80f8608889
# ╠═d5606f38-9231-463b-93fa-b74a7651fbcf
# ╠═a76601ad-6e79-4b07-9e99-5285a5e090db
# ╠═6dcd45b7-fdd4-4e7e-8ec6-e1efae25c78d
# ╠═b201e2fe-e14a-4d76-8266-181a40b3c89a
# ╠═9a9d6b70-8914-4aa1-bfde-57620f25ce23
# ╠═078571af-5e5b-4ba2-86b5-50949f16b737
# ╟─eefd1fc3-c0e3-4cfb-b42b-fb2ead981b33
# ╟─10a9dbab-eb69-4489-9f50-18a915468a78
# ╟─ce7de978-7e46-4851-a94f-8f0be028a438
# ╟─e675f460-d681-40b4-8479-2c722d7a9b5f
# ╟─5e4a93ff-fe0a-4dda-a979-ee6add0ba3d7
# ╟─6a618dfa-5123-47bf-9ac2-5d7fc5c32fe1
# ╟─9df0a409-1c37-4ea9-b935-f0b19ef842d9
# ╟─a9d23ab0-216b-4d5e-82f0-7933d5dcb552
# ╟─d0cef261-f6f0-4692-b217-e7e7ce29644f
# ╟─e8955866-0d1c-4423-89c6-cb2c521b3367
# ╟─0c551855-1ca7-4713-a916-cf97f590157d
# ╟─db5d1895-19e5-4457-b999-37a10c4fb708
# ╟─75e78f81-8fe4-4c85-b335-029a28176f79
# ╟─568df3d4-16a2-42d7-bae3-c5ef69a0d967
# ╟─d40fb873-15c1-4ead-9bf7-9a9df9cb6190
# ╟─4a51e4b2-e388-4c8c-9e2a-b3d46fc6f651
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
