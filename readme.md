## LiBr (aqueous lithium bromide) Property Routines

 The property routines described here are for the aqueous salt mixture LiBr/H2O. These routines were developed primarily as a part of the Sorption Systems Consortium at the University of Maryland.
The basic rotuines are provided in an external library named LiBrSSC.dll. The routines available in this library are listed in Table 1. The units of the properties must be as specified.

Table 1. Summary of LiBrSSC Routines

|Functions|Units|Comments|
|-|-|-|
|Thermal conductivity|W/m-K|k = libr_k(x,T)|
|Dynamic viscosity|cP|μ = libr_μ(x,T)|
|Specific heat|J/g-K|cp = libr_cp(x,T,p) or cp_sat = libr_cp(x,T)|
|Enthalpy|J/g|h = libr_h(x,T,p) or h_sat = libr_h(x,T)|
|Chemical potential of water|J/g|μw = libr_uw(x,T, p)|
|Chemical potential of LiBr|J/g|μLiBr = libr_us(x,T, p)|
|Volume|m^3/kg|v = libr_v(x,T)|
|Entropy|J/g-K|s = libr_s(x,T,p) or s_sat = libr_s(x,T)|
|Saturation pressure|kPa|Psat = libr_p(x,T)|
|Saturation temperature|K|Tsat = libr_t(x,p)|
|Saturation mass fraction|-|Xsat = libr_x(T,p)|
|Crystallization temperature|°C|Tcryst = libr_tCryst(x)|
|Index of refraction|-|N = libr_refindex(x,T)|
|Partial mass Gibbs function|J/g|g,dgdx,mu_w,mu_s=libr_part_g(x, T, P )|
|Partial mass enthalpy|J/g|h,dhdx,h_w,h_s=libr_part_h(x, T, P )|
|Partial mass entropy|J/g-K|s,dsdx,s_w,s_s=libr_part_s(x, T, P )|
|Partial mass volume|$m^3/kg$|v,dvdx,v_w,v_s=libr_part_v(x, T, P )|
|Flashing process|-|q,T,x,hl,hv=libr_flash(x,h,p)|
|help function|||
|Temperature according Enthalpy|K|T=libr_xht(x,h,p)|
|Mass faction accroding Enthalpy|%|x=libr_htx(h,t,p)|
|Temperature according entropy|K|T=libr_xst(x,s,p)|
|Mass faction according entropy|%|x=libr_stx(s,t,p)|

See details on units in procedure description
Input Units: Tc - K, X – mass fraction of LiBr %, P – kPa


### Transport Properties
The transport property correlations included were derived from an examination of available data. This was done as
part of the SSC efforts but was unpublished [1]. Details of the correlations are summarized here

### Viscosity
$$lnμ=A_0+A_1X^2+\frac{B_0}{T}+\frac{B_1X^2}{T}+\frac{C_0}{T}+\frac{C_1X^2}{T}$$

where the units are μ (cP), T (K), X (mass fraction LiBr) and the constants are

|Coef0|Coef1|
|-|-|
|Ao = -2.3212641667148|A1 = 3.190587778753|
|Bo = -609.44957160372|B1 = 963.16370163469|
|Co = 372994.85578423|C1 = -35211.99698739|

This equation yields a $R^2$ value of 0.984073 across all of the data sets used. Example calculation: T = 25°C (298.15
K), X = 50% ➡ μ = 3.807 cP.

### Thermal Conductivity
$$k=A_0+A_1X+B_0T+B_1TX+C_0T^2+C_1T^2X+D_0T^3+D_1T^3X$$

where k = W/m-K, T = K & X = mass fraction LiBr

|Coef0|Coef1|
|-|-|
|Ao= -0.880453887702949 |A1 = 0.883985046484968|
|Bo = 0.00898659269884302 |B1 = -0.007666522227789178|
|Co = -1.55427759660091E-05| C1 = 1.38873506415764E-05|
|Do = 7.3203107999836E-09| D1 = -6.31953452062666E-09|

Our “best fit” produced an R2 value of 0.9844 and the most well behaved curve over a broad range of the
independent parameters. There were equations that had higher R2 values but the one shown above had the best
combination of properties. Example calculation: T = 25°C (298.15 K), X = 50%  k = 0.444 W/m-K.

### Index of Refraction
The index of refraction of aqueous lithium bromide can be calculated from the correlation of Bostick et al. [2]

$$N=N_1X^2+N_2X+N_3T+N_4$$

where


|i|$N_i$|
|-|-|
|1| 0.0000241|
|2| 0.00108|
|3| -0.000106|
|4| 1.3348|

Units: T - °C, x – mass fraction LiBr in %

### Crystallization Temperature

The crystallization temperature is calculated from a curve fit to the data of Boryta [3].
$$x=A_0+A_1T+A_2T^2$$

<table>
<tr><td>i</td><td colspan="3" >Ai</td></tr>
<tr><td></td><td>65.05 < x < 71.91</td><td>57.08 < x < 65.05</td><td>48.47 < x < 57.08</td></tr>
<tr><td>0</td><td>62.63716</td><td>56.95202</td><td>56.55952</td></tr>
<tr><td>1</td><td>0.04810823</td><td>0.05205944</td><td>0.2337275</td></tr>
<tr><td>2</td><td>0.00024301</td><td>0.00346278</td><td>0.00141297</td></tr>
</table>

Units: T - °C, x – mass fraction LiBr in %

### Thermodynamic Properties

The thermodynamic properties include volume, enthalpy, specific heat, entropy and chemical potential as well as the
saturation properties temperature and pressure. All of these properties are derived from a Gibbs function fitted to a
broad set of data for the mixture. The advantage of this approach is that any thermodynamic property of interest can
be easily derived from the Gibbs function. And a corollary is that all of the derived properties will be thermodynamically consistent (up to the precision of the calculation). The details of the Gibbs function are available
in the literature [4, 5]. It provides excellent fidelity with the available data over the full concentration range from
pure water up to crystallization and from 0 – 300°C.

### Enthalpy of Aqueous Lithium Bromide

<center>h = libr_h(T,x) or h = libr_h(T,x,P)</center>
In general, we expect three independent properties to be required to specify state for a binary mixture such as
aqueous lithium bromide. However, because the pressure effect is usually small, it is often omitted from property
correlations for this mixture. The SSC routines retain and calculate pressure effects. However, to simplify routine
use where the pressure effects can be ignored, the SSC enthalpy routine has an overload feature which calculates
(and uses) the saturation pressure for cases where the user does not input a pressure value. Example calculation: T
= 25°C (298.15K), X = 50% ➡ h = 52.92 J/g.

### Entropy of Aqueous Lithium Bromide
<center>s = libr_s(x,T,p) or s = libr_s(x,T,p)</center>
Example calculation: T = 25°C (298.15K), X = 50% ➡ s = 0.1853 J/g-K.

### Specific Heat of Aqueous Lithium Bromide
<center>cp = libr_cp(x,T) or cp = libr_cp(x,T,p)</center>
Example calculation: T = 25°C (298.15K), X = 50% ➡ cp = 2.151 J/g-K.

### Specific Volume of Aqueous Lithium Bromide
<center>v = libr_v(x,T)</center>
Example calculation: T = 25°C (298.15K), X = 50% ➡ v = 0.0006523 m3/kg

### Chemical Potential in Aqueous Lithium Bromide
<center> g,dgdx,mu_w,mu_s=libr_part_g(x,T,p)</center>

$$g_w=u_w=g-x(\frac{∂g}{∂x})_{p,T}$$
$$g_s=u_s=g+(100-x)\frac{∂g}{∂x}_{p,T}$$

Example calculation: T = 25°C (298.15K), X = 50% LiBr, P = 0.8071 kPa ➡ g = -2.337 J/g, dg/dx=3.785 J/g, gw = -
191.6 J/g, gs = 186.9 J/g.

One of the aspects of partial properties is that you can sum them up to obtain the mixture property as
$$g=\frac{(100-x)g_w+xg_s}{100}$$

### Partial Enthalpy in Aqueous Lithium Bromide
<center> h,dhdx,h_w,h_s=libr_part_h(x,T,p)</center>

$$h_w=h-x(\frac{∂h}{∂x})_{p,T}$$
$$h_s=h+(100-x)\frac{∂h}{∂x}_{p,T}$$

Example calculation: T = 25°C (298.15K), X = 50% LiBr, P = 0.8071 kPa ➡ h = 52.92 J/g, dh/dx=1.944 J/g, hw = -
44.25 J/g, hs = 150.1 J/g.

One of the aspects of partial properties is that you can sum them up to obtain the mixture property as
$$h=\frac{(100-x)h_w+xh_s}{100}$$

### Partial Entropy in Aqueous Lithium Bromide
<center> s,dsdx,s_w,s_s=libr_part_s(x,T,p)</center>

$$s_w=s-x(\frac{∂s}{∂x})_{p,T}$$
$$s_s=s+(100-x)\frac{∂s}{∂x}_{p,T}$$

Example calculation: T = 25°C (298.15K), X = 50 % LiBr, P = 0.8071 kPa ➡ s = 0.1853 J/g-K, ds/dx=-0.006176
J/g-K, hw = 0.4942 J/g-K, ss = -0.1235 J/g-K.

One of the aspects of partial properties is that you can sum them up to obtain the mixture property as
$$s=\frac{(100-x)s_w+xs_s}{100}$$

### Partial Volume in Aqueous Lithium Bromide
<center> v,dvdx,v_w,v_s=libr_part_v(x,T)</center>

$$v_w=v-x(\frac{∂v}{∂x})_{p,T}$$
$$v_s=v+(100-x)\frac{∂v}{∂x}_{p,T}$$

Example calculation: T = 25°C (298.15K), X = 50% LiBr ➡ v = 0.6523 cm3/g, dv/dx=-0.006976 cm3/g, vw = 1.001
cm3/g, vs = 0.3033 cm3/g.

One of the aspects of partial properties is that you can sum them up to obtain the mixture property as
$$v=\frac{(100-x)v_w+xv_s}{100}$$

### Saturation Properties of Aqueous Lithium Bromide

* Psat = libr_p(x,T)
* Tsat = libr_t(x,p)
* Xsat = libr_x(T,p)

Example calculation: T = 25°C, X = 50% ➡ P = 0.8052 kPa

### Flashing of Aqueous Lithium bromide
<center>q,T,x,hl,hv=LibrFlash(x,h,p)</center>
This is a utility routine to automate the mass and energy balances involved in flashing through a valve or other
pressure restriction. It assumes an adiabatic process. Inputs are enthalpy (h), LiBr mass fraction and the pressure at
the outlet of the restriction. Outputs are the vapor quality (q), temperature (T), liquid mass fraction (x), liquid
enthalpy (hl) and vapor enthalpy (hv).
Example calculation: h = 74.54 J/g, P = 0.7 kPa, x_in = 50%  ➡ q = 0.009845, T = 23.33°C, xl = 50.49%, hl = 50.38 J/g,
hv = 2504.34 J/g.

Note: The overall mass, LiBr mass and energy balances are all satisfied

### References
1. SSC, Transport Property Data for Aqueous Lithium Bromide, in SSC
Unpublished Report1998.
2. Bostick, D.A., Klatt, L.N., Perez-Blanco, H., Fiber optics
refractometer for absorption machines, 1987.
3. Boryta, D.A., Solubility of Lithium Bromide in Water Between -50°C and
100°C (45 to 70% Lithium Bromide). J. Chem. Eng. Data,, 1970. 15(1): p.
142-144.
4. Yuan, Z. and K.E. Herold, Thermodynamic properties of aqueous lithium
bromide using a multiproperty free energy correlation. HVAC&R Research,
2005 11(3): p. 377-393.
5. Yuan, Z. and K.E. Herold, Specific heat measurements on aqueous lithium
bromide. HVAC&R Research, 2005. 11(3): p. 361-375.
