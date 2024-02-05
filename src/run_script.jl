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

Si_QE_input_path = joinpath(@__DIR__, "test.in")
Si_psp_path = joinpath(@__DIR__, artifact"pd_nc_sr_lda_standard_0.4.1_upf/Si.upf")

kgrid = [6, 6, 6]
Ecut = 60.0u"Ry"

""" Build a DFTK computation (system and calculator) from a 
QuantumEspresso input file.
"""
function build_computation_from_QE(QE_input_path, psp_path; temperature=1e-4)
	system = load_system(QE_input_path)
	system = attach_psp(system; Dict(symbol => psp_path))
	
	if functional == "lda"
		model_kwargs = (; functionals = [:lda_x, :lda_c_pw], temperature)
	elseif functional == "gga"
		model_kwargs = (; functionals = [:gga_x_pbe, :gga_c_pbe], temperature)
	else
		throw(ErrorException("Unrecognized class of functionals."))
	end
	basis_kwargs = (; kgrid, Ecut)
	scf_kwargs = (; is_converged = DFTK.ScfConvergenceEnergy(austrip(1e-9u"Ry")))
	calculator = DFTKCalculator(; model_kwargs, basis_kwargs, scf_kwargs, verbose=true)
	(; system, calculator)
end
		
""" Run energy computation for using DFTK. 
"""
function run_DFTK_energy(system, calculator)
	energy = AtomsCalculators.potential_energy(system, calculator)
	(; calculator.state)
end

function run_DFTK_vc_relax(system, calculator)
        linesearch =  BackTracking(c_1= 1e-4, ρ_hi= 0.8, ρ_lo= 0.1, order=2, maxstep=Inf)
        solver = OptimizationOptimJL.LBFGS(; linesearch)
        optim_options = (; solver, f_tol=1e-10, g_tol=1e-5, iterations=30,
                         show_trace=true, store_trace = true, allow_f_increases=true)
        results = minimize_energy!(system, calculator; procedure = "vc_relax", solver, optim_options...)
        (; results, calculator.state)
end

function run_QE(QE_input_path, psp_path)
	psp_folder, psp_filename = splitdir(psp_path)
	# Inject correct values for psp and paths.
	sed_regex = """s@\\(.*pseudo_dir = \\)\\(.*\\)\$@\\1"$(psp_folder)"/@g"""
	input_stream = read(pipeline(
		`sed -e "$sed_regex" $QE_input_path`,
		`perl -0777 -pe "s/.*(ATOMIC_SPECIES.*\n.[^0-9]*[0-9]*\.[0-9]*)(.*)/\1 $psp_filename/g"`), String)
	read(pipeline(`echo $input_stream`, `pw.x`), String)
end
