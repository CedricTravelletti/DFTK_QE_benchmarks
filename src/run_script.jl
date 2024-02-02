using DFTK
using AtomsCalculators
using Unitful
using UnitfulAtomic
using OptimizationOptimJL
using LazyArtifacts
using AtomsIO
using AtomsIOPython

using GeometryOptimization
using LineSearches

Si_QE_input_path = joinpath(@__DIR__, "si_relax.in")
Si_psp_path = joinpath(@__DIR__, artifact"pd_nc_sr_lda_standard_0.4.1_upf/Si.upf")


""" Run energy computation for Si using DFTK. 
"""
function run_Si_DFTK_energy(QE_input_path, psp_path)
	# Load system from QuantumEspresso file.
	system = load_system(QE_input_path)
	# PseudoDojo psp.
	system = attach_psp(system; Si=psp_path)
	
	model_kwargs = (; functionals = [:gga_x_pbe, :gga_c_pbe], temperature=1e-4)
	basis_kwargs = (; kgrid = [6, 6, 6], Ecut = 60.0u"Ry")
	scf_kwargs = (; is_converged = DFTK.ScfConvergenceEnergy(austrip(1e-9u"Ry")))
	calculator = DFTKCalculator(; model_kwargs, basis_kwargs, scf_kwargs, verbose=true)
	
	energy = AtomsCalculators.potential_energy(system, calculator)
	(; calculator.state)
end

function run_Zn_DFTK_vc_relax(QE_input_path, psp_path)
	system = load_system(QE_input_path)
	system = attach_psp(system; Si=psp_path)
        
        model_kwargs = (; functionals = [:lda_x, :lda_c_pw], temperature=1e-4)
        basis_kwargs = (; kgrid = [6, 6, 4], Ecut = 40.0u"Ry")
        scf_kwargs = (; is_converged = DFTK.ScfConvergenceEnergy(austrip(1e-9u"Ry")))
        calculator = DFTKCalculator(; model_kwargs, basis_kwargs, scf_kwargs, verbose=true)

        linesearch =  BackTracking(c_1= 1e-4, ρ_hi= 0.8, ρ_lo= 0.1, order=2, maxstep=Inf)
        solver = OptimizationOptimJL.LBFGS(; linesearch)
        optim_options = (; solver, f_tol=1e-10, g_tol=1e-5, iterations=30,
                         show_trace=true, store_trace = true, allow_f_increases=true)
        results = minimize_energy!(system, calculator; procedure = "vc_relax", solver, optim_options...)
        (; results, calculator.state)

function run_Si_QE_energy(Si_QE_input_path, Si_psp_path)
	pseudo_folder, pseudo_filename = splitdir(Si_psp_path)
	# Inject the pseudopotentials in the QE input file and run.
	input_stream = read(pipeline(
			`sed "s@%PSEUDO_DIR_PLACEHOLDER@\"$pseudo_folder/\"@g" $Si_QE_input_path`,
			`sed "s@%PSEUDO_FILENAME_PLACEHOLDER@$pseudo_filename@g"`), String)
	read(pipeline(`echo $input_stream`, `pw.x`), String)
end

