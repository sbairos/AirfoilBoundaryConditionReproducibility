using Trixi
using Trixi2Vtk
using OrdinaryDiffEq
using Plots

###############################################################################
# problem statment for compressible Euler equations of gas dynamics

equations = CompressibleEulerEquations2D(1.4)

# initial conditions for density, velocity, pressure
@inline function bc_mach2_flow(x, t, equations::CompressibleEulerEquations2D)
    rho_freestream = 1.4 # set the freestream flow parameters
    v1 = 80.0
    v2 = 0.0
    p_freestream = 1.0

    prim = SVector(rho_freestream, v1, v2, p_freestream)
    return prim2cons(prim, equations) # convert primitive to conservative variables
end
initial_condition = bc_mach2_flow

@inline function zero_dirichlet_bc(x, t, equations::CompressibleEulerEquations2D)
    prim = SVector(1.4, 0.0, 0.0, 1.0)
    return prim2cons(prim, equations)
end

@inline function bc_supersonic_inflow(u_inner, normal_direction::AbstractVector, x, t, surface_flux_function,
                                    equations::CompressibleEulerEquations2D)
    # Supersonic inflow boundary condition
    u_boundary = bc_mach2_flow(x, t, equations)
    flux = Trixi.flux(u_boundary, normal_direction, equations)
    return flux
end

@inline function bc_supersonic_outflow(u_inner, normal_direction::AbstractVector, x, t, surface_flux_function,
                                    equations::CompressibleEulerEquations2D)
    # Supersonic outflow boundary condition
    flux = Trixi.flux(u_inner, normal_direction, equations)
    return flux
end

###############################################################################
#  spatial discretization  of the equations using DG

polydeg = 3
surface_flux = flux_lax_friedrichs
volume_flux = flux_ranocha
basis = LobattoLegendreBasis(polydeg)
shock_indicator = IndicatorHennemannGassner(equations, basis, alpha_max = 0.5, alpha_min = 0.001, 
                                            alpha_smooth = true, variable = density_pressure)

volume_integral = VolumeIntegralShockCapturingHG(shock_indicator;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)
                                                
# Solver: Discontinuous Galerkin Spectral Element Method
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux, volume_integral = volume_integral)

###############################################################################
#  mesh from a gmsh geometry and rhs for semi-discretization 

meshSize = 0.1
geo_file = "NACA6412airfoil.geo"
mesh_file = string(geo_file,".inp")
run(`gmsh -setnumber meshSize $meshSize -2 -o $mesh_file $geo_file`)

mesh = P4estMesh{2}(mesh_file, polydeg = polydeg, 
                    boundary_symbols = [:PhysicalLine1, :PhysicalLine2, :PhysicalLine3, :PhysicalLine4])

boundary_conditions = Dict(:PhysicalLine1 => bc_supersonic_inflow, # left boundary
                    :PhysicalLine2 => bc_supersonic_outflow, # right boundary
                    :PhysicalLine3 => BoundaryConditionDirichlet(zero_dirichlet_bc), # Top and bottom boundary 
                    :PhysicalLine4 => BoundaryConditionDirichlet(zero_dirichlet_bc)) # Airfoil

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, boundary_conditions = boundary_conditions)

tspan = (0.0, 0.6) # time span to simulate
ode = semidiscretize(semi, tspan)

###############################################################################
#  callback functions during the time integration

summary_callback = SummaryCallback()

analysis_callback = AnalysisCallback(semi, interval = 100 )

stepsize_callback = StepsizeCallback(cfl = 0.8)

save_callback = SaveSolutionCallback(dt = 0.006, save_initial_solution=true, save_final_solution=true)

###############################################################################
# Run the simulation based on a Runge-Kutta method

callbacks = CallbackSet(summary_callback, analysis_callback, stepsize_callback, save_callback )

# positivity limiter necessary for this example with strong shocks. 
stage_limiter! = PositivityPreservingLimiterZhangShu(thresholds = (5.0e-7, 1.0e-6),
                                                    variables = (pressure, Trixi.density))

sol = solve(ode, SSPRK33(; thread = OrdinaryDiffEq.True(), stage_limiter! = stage_limiter!);
            dt = 1.0, save_everystep = false, callback = callbacks);

summary_callback() # print the timer summary

###############################################################################
# Basic plots

pd = PlotData2D(sol)
plot(pd, size=[1400, 800])
plot!(getmesh(pd))

trixi2vtk("out/solution*.h5", output_directory="out/")