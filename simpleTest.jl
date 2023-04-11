using JutulDarcy, Jutul, MultiComponentFlash
g = CartesianMesh((3, 1), (3.0, 1.0)) # to make porevolumes 1
geo = tpfv_geometry(g)
poro = [1.0, 0.1, 1.0]
G = discretized_domain_tpfv_flow(geo, porosity = poro)

nc = number_of_cells(G)

timesteps = [1.0, 2.0] .* 3600*24
inj = 1
prod = nc
#G.grid.pore_volumes[inj] *= 1000
#G.grid.pore_volumes[prod] *= 1000
co2 = MolecularProperty(0.0440, 7.38e6, 304.1, 9.412e-5, 0.224)
c1 = MolecularProperty(0.0160, 4.60e6, 190.6, 9.863e-5, 0.011)
c10 = MolecularProperty(0.0142, 2.10e6, 617.7, 6.098e-4, 0.488)

z0 = [0.5, 0.3, 0.2]
zi = [0.99, 0.01-1e-3, 1e-3]
mixture = MultiComponentMixture([co2, c1, c10], names = ["CO2", "C1", "C10"])

p0 = 75e5
T0 = 423.25

n = length(z0)
eos = GenericCubicEOS(mixture)
nc = number_of_cells(G)
L, V = LiquidPhase(), VaporPhase()
# Define system and realize on grid
sys = MultiPhaseCompositionalSystemLV(eos, (L, V))
model = SimulationModel(G, sys)

push!(model.output_variables, :Saturations)
push!(model.output_variables, :PhaseMassDensities)
push!(model.output_variables, :PhaseViscosities)

kr = BrooksCoreyRelPerm(sys, [2, 3])
s = model.secondary_variables
s[:RelativePermeabilities] = kr
parameters = setup_parameters(model, Temperature = T0)

tot_time = sum(timesteps)
forces = setup_forces(model)

p_init = repeat([p0], nc)
p_init[inj] = 2*p0
p_init[prod] = p0/2

z_init = repeat(z0, 1, nc)
z_init[:, inj] .= zi

# State is dict with pressure in each cell
init = Dict(:Pressure => p_init, :OverallMoleFractions => z_init)
if isnothing(parameters)
    parameters = setup_parameters(model)
end
state0 = setup_state(model, init)

states, = simulate(state0, model, timesteps, forces = forces, parameters = parameters);

## more information
sim = Simulator(model, state0 = state0)
simulate!(sim, timesteps, forces = forces, parameters = parameters, max_timestep_cuts = 0, max_nonlinear_iterations = 0);

# block-structure