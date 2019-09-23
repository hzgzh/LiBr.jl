"""
# libr_
Li-Br solution properties for absorb chiller
See also: Thermodynamic Properties of Aqueous Lithium Bromide Using a Multiproperty Free Energy Correlation Zhe Yuan Keith E. Herold Member ASHRAE Received June 10, 2004; accepted January 7, 2005 A

## Export function


"""
module LiBr

export libr_k,
       libr_h,
       libr_s,
       libr_cp,
       libr_μ,
       libr_uw,
       libr_us,
       libr_v,
       libr_s,
       libr_p,
       libr_t,
       libr_x,
       libr_tcryst,
       libr_flash,
       libr_refindex,
       libr_part_g,
       libr_part_h,
       libr_part_s,
       libr_part_v


using Roots
using XSteam
using ForwardDiff
import Base.sum

A = [5.506219979e3, 5.213228937e2, 7.774930356, -4.575233382e-2, -5.792935726e2]
B = [
    1.452749674e2,
    -4.984840771e-1,
    8.83691918e-2,
    -4.870995781e-4,
    -2.905161205,
]
C = [
    2.648364473E-2,
    -2.311041091E-3,
    7.559736620E-6,
    -3.763934193E-8,
    1.176240649E-3,
]
D = [-8.526516950e-6, 1.320154794e-6, 2.791995438e-11, 0, -8.511514931e-7]
E = [-3.840447174e-11, 2.625469387e-11]
F = [-5.159906276e1, 1.114573398]
L = [
    -2.183429482e3,
    -1.266985094e2,
    -2.364551372,
    1.389414858e-2,
    1.583405426e2,
]
M = [
    -2.267095847e1,
    2.983764494e-1,
    -1.259393234e-2,
    6.849632068e-5,
    2.767986853e-1,
]
V = [
    1.176741611e-3,
    -1.002511661e-5,
    -1.695735875e-8,
    -1.497186905e-6,
    2.538176345e-8,
    5.815811591e-11,
    3.057997846e-9,
    -5.129589007e-11,
]
ce = [0.0, 1.0, 2.0, 3.0, 1.1]
c0 = [1.0, 2.0, 3.0, 1.1]
c1 = [0.0, 1.0, 2.0, 0.1]

T0 = 220.0
sum(A, x, xc, y = 1.0, yc = [1.0]) = Base.sum(@. A * x^xc * y^yc)

function g(x, T, p = 1.0)

    re = [
        sum(A, x, ce),
        T * sum(B, x, ce),
        T^2.0 * sum(C, x, ce),
        T^3.0 * sum(D, x, ce),
        T^4.0 * sum(E, x, [0.0, 1.0]),
        sum(F, x, [0.0, 1.0]) / (T - T0),
        p *
        sum(
            V,
            x,
            [0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0],
            T,
            [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0],
        ),
        log(T) * sum(L, x, ce),
        T * log(T) * sum(M, x, [0.0, 1.0, 2.0, 3.0, 1.1]),
    ]
    sum(re)
end
∂g∂x(x, T, p = 1.0) = ForwardDiff.derivative(x -> g(x, T, p), x)
∂²g∂x²(x, T, p = 1.0) = ForwardDiff.derivative(x -> ∂g∂x(x, T, p), x)
∂g∂t(x, T, p = 1.0) = ForwardDiff.derivative(T -> g(x, T, p), T)
∂g∂p(x, T, p) = ForwardDiff.derivative(p -> g(x, T, p), p)
∂²g∂t²(x, T, p) = ForwardDiff.derivative(T -> ∂g∂t(x, T, p), T)
∂h∂x(x, T, p) = ForwardDiff.derivative(x -> libr_h(x, T, p), x)

libr_s(x, T, p = 1.0) = -∂g∂t(x, T, p)
libr_v(x, T, p = 1.0) = ∂g∂p(x, T, p)
libr_cp(x, T, p) = -T * ∂²g∂t²(x, T, p)
libr_h(x, T, p) = g(x, T, p) + T * libr_s(x, T, p)
libr_uw(x, T, p = 1.0) = g(x, T, p) - x * ∂g∂x(x, T, p)
libr_us(x, T, p) = g(x, T, p) + (100 - x) * ∂g∂x(x, T, p)
libr_u(x, T, p) = libr_h(x, T, p) - p * libr_v(x, T)
hw(x, T, p) = libr_h(x, T, p) - x * ∂h∂x(x, T, p)
hs(x, T, p) = libr_h(x, T, p) + (100.0 - x) * ∂h∂x(x, T, p)
qd(x, T, p) = libr_h(0.0, T, p) - hw(x, T, p)

function libr_refindex(x, T)

    N = [0.0000241, 0.00108, -0.000106, 1.3348]
    N[1] * x^2 + N[2] * x + N[3] * T + N[4]
end

function libr_tcryst(x)

    a0 = [62.63716, 0.04810823, 0.00024301]
    a1 = [56.95202, 0.05205944, 0.00346278]
    a2 = [56.55952, 0.2337275, 0.00141297]
    if x > 48.47 && x <= 57.08
        return (-a2[2] + sqrt(a2[2] * a2[2] - 4.0 * a2[3] * (a2[1] - x))) /
               2.0 / a2[3]
    elseif x > 57.08 && x <= 65.05
        return (-a1[2] + sqrt(a1[2] * a1[2] - 4.0 * a1[3] * (a1[1] - x))) /
               2.0 / a1[3]
    elseif x > 65.05 && x <= 71.91
        return (-a0[2] + sqrt(a0[2] * a0[2] - 4.0 * a0[3] * (a0[1] - x))) /
               2.0 / a0[3]
    else
        throw(DomainError("x should in  range of 48.47%-71.91%"))
    end
end

function libr_μ(x, T)

    x = x / 100
    a = [-2.3212641667148, 3.190587778753]
    b = [-609.44957160372, 963.16370163469]
    c = [372994.85578423, -35211.99698739]
    d = a[1] + a[2] * x^2.0 + b[1] / T + b[2] * x^2.0 / T + c[1] / T^2.0 +
        c[2] * x^2.0 / T^2.0
    return exp(d)
end

function libr_k(x, T)

    x = x / 100.0
    a = [-0.880453887702949, 0.883985046484968]
    b = [0.00898659269884302, -0.007666522227789178]
    c = [-1.55427759660091e-5, 1.38873506415764e-5]
    d = [7.3203107999836e-9, -6.31953452062666e-9]
    return a[1] + a[2] * x + b[1] * T + b[2] * T * x + c[1] * T * T +
           c[2] * T * T * x + d[1] * T^3.0 + d[2] * T^3.0 * x
end

fun(x, T, p) =
    h_pT(p / 100.0, T - 273.15) - (T) * s_pT(p / 100.0, T - 273.15) -
    libr_uw(x, T, p)

function libr_p(x, T)
    f(p) = fun(x, T, p)
    fzero(f, 0.8)
end

function libr_x(T, p)
    f(x) = fun(x, T, p)
    fzero(f, 50.0)
end

function libr_t(x, p)
    f(T) = fun(x, T, p)
    fzero(f, 330.0)
end

function flashfun(xl, x, h1, p)
    hv = hV_p(p / 100)
    T = libr_t(xl, p)
    h2 = libr_h(xl, T, p)
    hv * (1.0 - x / xl) + x / xl * h2 - h1
end

function libr_flash(x, h, p)
    if p < 0.613
        throw(DomainError(" error: p<0.613kPa"))
    end
    f(xt) = flashfun(xt, x, h, p)
    xt = fzero(f, x)
    hv = hV_p(p / 100)
    q = 1.0 - x / xt
    T2 = libr_t(xt, p)
    hl = (h - q * hv) / (1 - q)
    (q, T2, xt, hl, hv)
end

function libr_part_g(x, T, p)
    gx = g(x, T, p)
    g_x = ∂g∂x(x, T, p)
    u_w = libr_uw(x, T, p)
    u_s = libr_us(x, T, p)
    gx, g_x, u_w, u_s
end

function libr_part_h(x, T, p)
    h = libr_h(x, T, p)
    ∂h∂x = ForwardDiff.derivative(x -> libr_h(x, T, p), Float64(x))
    hw = h - x * ∂h∂x
    hs = h + (100 - x) * ∂h∂x
    h, ∂h∂x, hw, hs
end

function libr_part_s(x, T, p)
    s = libr_s(x, T, p)
    ∂s∂x = ForwardDiff.derivative(x -> libr_s(x, T, p), x)
    sw = s - x * ∂s∂x
    ss = s + (100 - x) * ∂s∂x
    s, ∂s∂x, sw, ss
end

function libr_part_v(x, T)
    v = libr_v(x, T)
    ∂v∂x = Calculus.derivative(x -> libr_v(x, T), x)
    vw = v - x * ∂v∂x
    vs = v + (100 - x) * ∂v∂x
    v, ∂v∂x, vw, vs
end


end #module
