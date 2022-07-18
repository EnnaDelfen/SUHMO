# Model equations

Go back to the [main documentation page](https://ennadelfen.github.io/SUHMO/index)

## Introduction 
SUHMO stands for SUblglacial Hydrology MOdel, and is an Adaptive Mesh Refinement (AMR) model based on the Chombo[^1] software framework. We extend the model proposed by Sommers et al. (2018)[^2] with a few changes to accommodate the transition from unresolved to resolved flow features. We handle the strong nonlinearities present in the equations by resorting to an efficient nonlinear Full Approximation Scheme multigrid (FAS-MG) algorithm. We outline the details of the algorithm next.

## Model Equations
The governing equation set starts with a two-dimensional expression for the conservation of mass – assuming we are dealing with an incompressible fluid:

$$\frac{\partial b}{\partial t} + \frac{\partial b_e}{\partial t} + \nabla \cdot \mathbf{q} = \frac{\dot{m}}{\rho_w} + e_s$$

where $b$ is the subglacial water gap height (m), $b_e$ is the volume of water stored englacially per unit area of bed (m), $\mathbf{q}$ is the gap-integrated basal water flux (m$^2$ s$^{-1}$), $\dot{m}$ is the melt rate (kg m$^{-2}$ s$^{-1}$), $\rho_w$ is the density of water (kg m$^{-3}$), and $e_s$ encompasses all external sources of meltwater (produced englacially or surface meltwater, for example) (m s$^{-1}$).

An approximate momentum equation for water velocity integrated over the gap height gives rise to an expression for the water flux, based on equations developed for flow in rock fractures[^3]:

$$\mathbf{q} = \frac{-b^3 g}{12 \nu (1 + \omega Re)} \nabla h$$

where $g$ is the gravitational acceleration (m s$^{-2}$), $\nu$ is the kinematic viscosity of water (m$^{2}$ s$^{-1}$), $\omega$ is a dimensionless parameter controlling the nonlinear transition from laminar to turbulent flow and $Re$ is the Reynolds number. The hydraulic head $h$ (m) is defined:

$$ h=\frac{P_w}{\rho_w g} + z_b$$

where $P_w$ is the water pressure (Pa) and $z_b$ is the bed elevation (m). The Reynolds number follows a classical definition:

$$ Re = \frac{|v|b}{\nu} = \frac{|\mathbf{q}|}{\nu} $$

where $v$ is the average flow velocity across the gap height. 

The melt rate $\dot{m}$ includes heat produced at the bed (geothermal flux and frictional heat due to sliding over the bed) along with heat generated through internal dissipation (mechanical energy converted to thermal energy by the flow), which is effectively melting the drainage system's walls and ceiling: 

$$ \dot{m} = \frac{1}{L} (G + |\mathbf{u}_b \cdot \mathbf{\tau}_b|
    - \rho_w g \mathbf{q} \cdot \nabla h
    + c_t c_w \rho_w \mathbf{q} \cdot \nabla P_w) $$
    
where $L$ is the latent heat of fusion of water (J kg$^{-1}$), $G$ is the geothermal flux (W m$^{-2}$), $\mathbf{u}_b$ is the ice basal velocity vector (m s$^{-1}$), $\mathbf{\tau}_b$ is the stress exerted by the bed onto the ice (Pa), $c_t$ is the change in pressure melting point with temperature (K Pa$^{-1}$), and $c_w$ is the heat capacity of water (J kg$^{-1}$ K$^{-1}$). We note that the last term takes into account the changes in sensible heat due to pressure melting point variations. This term is often considered negligible and dropped from similar models[^4].

Finally, the effective drainage-system capacity $b'$ evolves according to opening and closure terms that are typically model-specific. Opening can be due to melt and sliding over bumps on the bed, while closing is solely due to ice creep:

$$  \frac{\partial b'}{\partial t} = \frac{\dot{m}}{\rho_i} + \beta u_b - A |P_i-P_w|^{n-1} (P_i - P_w) l_c $$

where $\rho_i$ is the ice density (kg m$^3$), $u_b$ is the magnitude of the sliding velocity (m s$^{-1}$), $A$ is the ice-flow parameter (Pa$^{-3}$ s$^{-1}$), $n$ is the flow-law exponent (typically, $n=3$) and $P_i$ is the ice overburden pressure (Pa).
The parameter $\beta= max((b_r -b)/{l_r}, 0)$ is dimensionless; it governs opening by sliding and is a function of the typical bed bump height ($b_r$) and bump spacing ($l_r$) in such a way that opening by sliding only occurs where the gap height is less than the typical local bump height.
The quantity $l_c$ is the creep length scale, which is generally equal to $b$.

We do not allow for the drainage space to be partially filled, such that $b=b'$ always. These equations can then be combined to produce an equation for the evolution of the hydraulic head: 

$$ \nabla \cdot \Big[ \frac{-b^3g}{12 \nu (1 + \omega Re)} \nabla h \Big] + \frac{\partial b_e}{\partial t} = \dot{m} \Big[ \frac{1}{\rho_w} - \frac{1}{\rho_i} \Big] + A |P_i - P_w|^{n-1}(P_i - P_w) l_c -\beta u_b + e_s $$




[^1]:Adams, M., Colella, P., Graves, D. T., Johnson, J. N., Keen, N. D., Ligocki, T. J., Martin, D. F., McCorquodale, P. W., Modiano, D., Schwartz, P. O., Sternberg, T. D., and Van Straalen, B.: Chombo Software Package for AMR Applications - Design Document, Tech. Rep. LBNL-6616E, Lawrence Berkeley National Laboratory, [link](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations), 2001-2021.
[^2]:Sommers, A., Rajaram, H., and Morlighem, M.: SHAKTI: Subglacial Hydrology and Kinetic, Transient Interactions v1.0., Geoscientific Model Development, 11, 2955–2974, [doi](https://doi.org/10.5194/gmd-11-2955-2018), 2018.
[^3]:Zimmerman, R. W., Al-Yaarubi, A., Pain, C. C., and Grattoni, C. A.: Non-linear regimes of fluid flow in rock fractures, International Journal of Rock Mechanics and Mining Sciences, 41, 163–169, [doi](https://doi.org/10.1016/j.ijrmms.2003.12.045), 2004.
[^4]:de Fleurian, B., Werder, M. A., Beyer, S., Brinkerhoff, D. J., Delaney, I., Dow, C. F., Downs, J., Gagliardini, O., Hoffman, M. J., Hooke, R. L., et al.: SHMIP The subglacial hydrology model intercomparison Project, Journal of Glaciology, 64, 897–916, [doi](https://doi.org/10.1017/jog.2018.78), 2018.

