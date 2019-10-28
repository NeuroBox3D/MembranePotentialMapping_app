--------------------------------------------------------------
--  Example script for simulation on 1d full neuron cell    --
--  solving the HH cable problem with 5 synapses            --
--  The voltage at each (x,y,z) is a specified from a .dat  --
--  file that is read in                                    --
--                                                          --
--  Author(s): James Rosado,  Stephan Grein                 --
--                                                          --
--  Date: Sept-2019                                         --
--------------------------------------------------------------

-- ug specific loads 
SetOutputProfileStats(false)
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- require plugins to be available
AssertPluginsLoaded({"cable_neuron", "neuro_collection", "MembranePotentialMapping"})

-------------------------------------------------------------------------------
-- parameters
-------------------------------------------------------------------------------
dim = util.GetParamNumber("-dim", 3)
gridName1d = util.GetParam("-grid1d", "my_cell.ugx")
numRefs = util.GetParamNumber("-numRefs", 0)
gridSyn = string.sub(gridName1d, 1, string.len(gridName1d) - 4) .. "_syns.ugx" 
vdccMode = util.GetParamBool("-vdccMode", false) 
-- time parameters in [s]
dt1d = util.GetParamNumber("-dt1d", 1e-5)
dt3d = util.GetParamNumber("-dt3d", 1e-2)
dt3dStart = util.GetParamNumber("-dt3dstart", dt3d)
endTime = util.GetParamNumber("-endTime", 1.0) 
-- specify "-verbose" to output linear solver convergence
verbose1d = util.HasParamOption("-verbose1d")
verbose3d = util.HasParamOption("-verbose3d")
-- vtk output?
-- generateVTKoutput = util.HasParamOption("-vtk")
generateVTKoutput =util.GetParam("-vtk",true)
pstep = util.GetParamNumber("-pstep", dt3d, "plotting interval")

-- file handling
filename = util.GetParam("-outName", "hybrid_test")
filename = filename.."/"

-- init ug 
InitUG(dim, AlgebraType("CPU", 1));

-------------------------------------------------------------------------------
-- solver
-------------------------------------------------------------------------------
solverID = util.GetParam("-solver", "GMG")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG"] = 0;
validSolverIDs["GS"] = 0;
validSolverIDs["ILU"] = 0;
if not validSolverIDs[solverID] then error("Unknown solver ID: " .. solverID) end

-- choose length of time step at the beginning
-- if not timeStepStart = 2^(-n)*timeStep, take nearest lower number of that form
function log2(x)
  return math.log(x) / math.log(2)
end

startLv =  math.ceil(log2(dt3d / dt3dStart))
dt3dStartNew = dt3d / math.pow(2, startLv)
if (math.abs(dt3dStartNew - dt3dStart) / dt3dStart > 1e-5) then 
  print("dt3dStart argument ("..dt3dStart..") was not admissible; taking "..dt3dStartNew.." instead.")
end
dt3dStart = dt3dStartNew

-------------------------------------------------------------------------------
-- print setup
-------------------------------------------------------------------------------
print("Chosen parameters:")
print("    grid       = " .. gridName1d)
print("    numRefs    = " .. numRefs)
print("    dt1d       = " .. dt1d)
print("    dt3d       = " .. dt3d)
print("    dt3dStart  = " .. dt3dStart)
print("    endTime    = " .. endTime)
print("    pstep      = " .. pstep)
print("    ions       = " .. tostring(withIons))
print("    solver     = " .. solverID)
print("    verbose1d  = " .. tostring(verbose1d))
print("    verbose3d  = " .. tostring(verbose3d))
print("    vtk        = " .. tostring(generateVTKoutput))

-------------------------------------------------------------------------------
-- setup synapses
-------------------------------------------------------------------------------
-- firing pattern of the synapse
syn_onset = {0.0, 0.0, 0.0, 0.0, 0.0}
syn_tau = 4e-4
syn_gMax = 1.2e-9
syn_revPot = 0.0

function file_exists(name)
  local f = io.open(name,"r")
  if f ~= nil then io.close(f) return true else return false end
end

-- check if grid version with synapses already exists
-- if so, just use it, otherwise, create it
if not file_exists(gridSyn) then
  -- synapse distributor only works in serial mode
  if NumProcs() > 1 then
    print("Cannot use SynapseDistributor in parallel atm. Please create synapse geometry in serial.")
    exit()
  end
  
  -- use synapse distributor to place synapses on the grid
  synDistr = SynapseDistributor(gridName1d)

  syn1 = AlphaSynapsePair()
  syn1:set_id(0)
  syn1:set_onset(syn_onset[1])
  syn1:set_tau(syn_tau)                   -- default value
  syn1:set_gMax(syn_gMax)                 -- default value
  syn1:set_reversal_potential(syn_revPot) -- default value
  --synDistr:place_synapse_at_coords({-19.9609e-6, 43.9832e-6, 5.82852e-6}, syn1:pre_synapse(), syn1:post_synapse())
  synDistr:place_synapse_at_coords({-86.74e-6, 590.09e-6, -10.55e-6}, syn1:pre_synapse(), syn1:post_synapse())
  
  syn2 = AlphaSynapsePair()
  syn2:set_id(1)
  syn2:set_onset(syn_onset[2])
  syn2:set_tau(syn_tau)                   -- default value
  syn2:set_gMax(syn_gMax)                 -- default value
  syn2:set_reversal_potential(syn_revPot) -- default value
  --synDistr:place_synapse_at_coords({-18.5196e-6, 40.3637e-6, 5.8769e-6}, syn2:pre_synapse(), syn2:post_synapse())
  synDistr:place_synapse_at_coords({-528.3e-6, 89.84e-6, -93.67e-6}, syn2:pre_synapse(), syn2:post_synapse())
  
  syn3 = AlphaSynapsePair()
  syn3:set_id(2)
  syn3:set_onset(syn_onset[3])
  syn3:set_tau(syn_tau)                   -- default value
  syn3:set_gMax(syn_gMax)                 -- default value
  syn3:set_reversal_potential(syn_revPot) -- default value
  --synDistr:place_synapse_at_coords({-18.6295e-6, 35.8769e-6, 5.87629e-6}, syn3:pre_synapse(), syn3:post_synapse())
  synDistr:place_synapse_at_coords({-404.99e-6, -367.3e-6, -55.98e-6}, syn3:pre_synapse(), syn3:post_synapse())
  
  syn4 = AlphaSynapsePair()
  syn4:set_id(3)
  syn4:set_onset(syn_onset[4])
  syn4:set_tau(syn_tau)                   -- default value
  syn4:set_gMax(syn_gMax)                 -- default value
  syn4:set_reversal_potential(syn_revPot) -- default value
  --synDistr:place_synapse_at_coords({-19.53e-6, 32.1163e-6, 5.87559e-6}, syn4:pre_synapse(), syn4:post_synapse())
  synDistr:place_synapse_at_coords({-447.6e-6, -110.29e-6, 47.81e-6}, syn4:pre_synapse(), syn4:post_synapse())
  
  syn5 = AlphaSynapsePair()
  syn5:set_id(4)
  syn5:set_onset(syn_onset[5])
  syn5:set_tau(syn_tau)                   -- default value
  syn5:set_gMax(syn_gMax)                 -- default value
  syn5:set_reversal_potential(syn_revPot) -- default value
  --synDistr:place_synapse_at_coords({-20.7218e-6, 28.1453e-6, 5.93687e-6}, syn5:pre_synapse(), syn5:post_synapse())
  synDistr:place_synapse_at_coords({265.79e-6, 325.2e-6, 39.88e-6}, syn5:pre_synapse(), syn5:post_synapse())
  
  exportSuccess = synDistr:export_grid(gridSyn)
  if not exportSuccess then
      print("Synapse distributor failed to export its grid to '" .. gridSyn .. "'.")
      exit()
  end
end

gridName1d = gridSyn

-------------------------------------------------------------------------------
-- biological settings
-------------------------------------------------------------------------------
-- settings are according to T. Branco

-- membrane conductances (in units of S/m^2)
g_k_ax = 400.0  -- axon
g_k_so = 200.0  -- soma
g_k_de = 30   -- dendrite

g_na_ax = 3.0e4
g_na_so = 1.5e3
g_na_de = 40.0

g_l_ax = 200.0
g_l_so = 1.0
g_l_de = 1.0

-- specific capacitance (in units of F/m^2)
spec_cap = 1.0e-2

-- resistivity (in units of Ohm m)
spec_res = 1.5

-- reversal potentials (in units of V)
e_k  = -0.09
e_na = 0.06
e_ca = 0.14

-- equilibrium concentrations (in units of mM)
-- comment: these concentrations will not yield Nernst potentials
-- as given above; pumps will have to be introduced to achieve this
-- in the case where Nernst potentials are calculated from concentrations!
k_out  = 4.0
na_out = 150.0
ca_out = 1.5

k_in   = 140.0
na_in  = 10.0
ca_in  = 5e-5

-- equilibrium potential (in units of V)
v_eq = -0.065

-- diffusion coefficients (in units of m^2/s)
diff_k  = 1.0e-9
diff_na = 1.0e-9
diff_ca = 2.2e-10

-- temperature in units of deg Celsius
temp = 37.0

-------------------------------------------------------------------------------
-- create 1d domain and approx space --
-------------------------------------------------------------------------------
neededSubsets1d = {"soma", "dendrite", "axon"}
dom1d = util.CreateDomain(gridName1d, 0, neededSubsets1d)
scale_domain(dom1d, 1e-6)

approxSpace1d = ApproximationSpace(dom1d)
approxSpace1d:add_fct("v", "Lagrange", 1)

approxSpace1d:init_levels();
approxSpace1d:init_surfaces();
approxSpace1d:init_top_surface();
approxSpace1d:print_layout_statistic()
approxSpace1d:print_statistic()

OrderCuthillMcKee(approxSpace1d, true);

-------------------------------------------------------------------------------
-- create 1d disc
-------------------------------------------------------------------------------
allSubsets = "soma, dendrite, axon"

-- cable equation
CE = CableEquation(allSubsets, false)
CE:set_spec_cap(spec_cap)
CE:set_spec_res(spec_res)
CE:set_rev_pot_k(e_k)
CE:set_rev_pot_na(e_na)
CE:set_rev_pot_ca(e_ca)
CE:set_k_out(k_out)
CE:set_na_out(na_out)
CE:set_ca_out(ca_out)
CE:set_diff_coeffs({diff_k, diff_na, diff_ca})
CE:set_temperature_celsius(temp)

-- Hodgkin and Huxley channels
HH = ChannelHHNernst("v, k, na", "axon")
HH = ChannelHH("v", allSubsets)
HH:set_conductances(g_k_ax, g_na_ax, "axon")
HH:set_conductances(g_k_so, g_na_so, "soma")
HH:set_conductances(g_k_de, g_na_de, "dendrite")

CE:add(HH)

-- leakage
tmp_fct = math.pow(2.3,(temp-23.0)/10.0)

leak = ChannelLeak("v", allSubsets)
leak:set_cond(g_l_ax*tmp_fct, "axon")
leak:set_rev_pot(-0.066148458, "axon")
leak:set_cond(g_l_so*tmp_fct, "soma")
leak:set_rev_pot(-0.030654022, "soma")
leak:set_cond(g_l_de*tmp_fct, "dendrite")
leak:set_rev_pot(-0.057803624, "dendrite")

CE:add(leak)

-- synapses
syn_handler = SynapseHandler()
syn_handler:set_ce_object(CE)
CE:set_synapse_handler(syn_handler)

-- domain discretization
domDisc1d = DomainDiscretization(approxSpace1d)
domDisc1d:add(CE)

assTuner = domDisc1d:ass_tuner()

-- setup time discretization --
timeDisc = ThetaTimeStep(domDisc1d)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler

-- create operator from discretization
linOp = AssembledLinearOperator(timeDisc)

-- Input Voltage Function
if vdccMode then
  function InputVm(x, y, z, t, si)
    -- Some random function to test this out
    return ((math.sin(x)*math.sin(y)*math.cos(z)) * si)-0.065
  end
  InputVmFunction = LuaUserNumber("InputVm")
else
  local mapper = Mapper()
  local voltageDataFile = "voltageData.csv"
  mapper:build_tree_from_file(voltageDataFile)
  function InputVm(x, y, z, t, si)
    return mapper:get_data_from_nn({x, y, z})
  end
  InputVmFunction = LuaUserNumber("InputVm")
end
---------------------------------------------------------------------------------

-------------------------------------------------------------------------------
-- solver setup --
-------------------------------------------------------------------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace1d)
dbgWriter:set_vtk_output(true)

-- linear solver --
convCheck = CompositeConvCheck(approxSpace1d, 1, 2e-26, 1e-08)
convCheck:set_component_check("v", 1e-21, 1e-12)
convCheck:set_verbose(verbose1d)

ilu = ILU()
solver = LinearSolver()
solver:set_preconditioner(ilu)
solver:set_convergence_check(convCheck)
--solver:set_debug(dbgWriter)

-------------------------------------------------------------------------------
-- solving
-------------------------------------------------------------------------------
-- get grid function
u = GridFunction(approxSpace1d)
b = GridFunction(approxSpace1d)

-- set initial value
InterpolateInner(InputVmFunction, u, "v", 0.0)

-- timestep in seconds
dt = dt3dStart
time = 0.0
step = 0

-- initial vtk output
if generateVTKoutput then
  out = VTKOutput()
  out:print(filename.."solution", u, step, time)
end

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

curr_dt = dt1d
dtred = 2

lv = 0
maxLv = 10
cb_counter = {}
cb_counter[lv] = 0
while endTime-time > 0.001*curr_dt do
  -- setup time Disc for old solutions and timestep
  timeDisc:prepare_step(solTimeSeries, curr_dt)
  
  -- reduce time step if cfl < curr_dt
  -- (this needs to be done AFTER prepare_step as channels are updated there)
  dtChanged = false
  cfl = CE:estimate_cfl_cond(solTimeSeries:latest())
  print("estimated CFL condition: dt < " .. cfl)
  while (curr_dt > cfl) do
    curr_dt = curr_dt/dtred
    
    if lv+1 > maxLv then
      print("Time step too small.")
      exit()
    end
    
    lv = lv + 1
    cb_counter[lv] = 0
    print("estimated CFL condition: dt < " .. cfl .. " - reducing time step to " .. curr_dt)
    dtChanged = true
  end
  
  -- increase time step if cfl > curr_dt / dtred (and if time is aligned with new bigger step size)
  while curr_dt*dtred < cfl and lv > 0 and cb_counter[lv] % (dtred) == 0 do
    curr_dt = curr_dt*dtred;
    lv = lv - 1
    cb_counter[lv] = cb_counter[lv] + cb_counter[lv+1]/dtred
    cb_counter[lv+1] = 0
    print ("estimated CFL condition: dt < " .. cfl .. " - increasing time step to " .. curr_dt)
    dtChanged = true
  end
  
  print("++++++ POINT IN TIME " .. math.floor((time+curr_dt)/curr_dt+0.5)*curr_dt .. " BEGIN ++++++")
  
  -- prepare again with new time step size
  if dtChanged == true then 
    timeDisc:prepare_step(solTimeSeries, curr_dt)
  end

  -- assemble linear problem
  matrixIsConst = time ~= 0.0 and dtChanged == false
  assTuner:set_matrix_is_const(matrixIsConst)
  AssembleLinearOperatorRhsAndSolution(linOp, u, b)
  
  -- apply linear solver
  ilu:set_disable_preprocessing(matrixIsConst)
  if ApplyLinearSolver(linOp, u, b, solver) == false then
    print("Could not apply linear solver.")
    exit()
  end
  
  -- update to new time
  time = solTimeSeries:time(0) + curr_dt
  
  -- vtk output
  if generateVTKoutput then
    if math.abs(time/pstep - math.floor(time/pstep+0.5)) < 1e-5 then 
      out:print(filename.."solution", u, math.floor(time/pstep+0.5), time)
    end
  end
  
  -- updte time series (reuse memory)
  oldestSol = solTimeSeries:oldest()
  VecScaleAssign(oldestSol, 1.0, u)
  solTimeSeries:push_discard_oldest(oldestSol, time)
  
  -- increment check-back counter
  cb_counter[lv] = cb_counter[lv] + 1
    
  print("++++++ POINT IN TIME  " .. math.floor(time/curr_dt+0.5)*curr_dt .. "s  END ++++++++");
end

-- end timeseries, produce gathering file
if (generateVTKoutput) then out:write_time_pvd(filename .. "solution", u) end
