using DFTK
using AtomsCalculators
using AtomsIO
using AtomsIOPython
using AtomsBase
using Unitful
using UnitfulAtomic


""" Build a DFTK computation (system and calculator) from a 
QuantumEspresso input file.
"""
function build_computation_from_QE(QE_input_path, psp_path;
				   atom_symbol, functional,
				   kgrid, Ecut, temperature=1e-4)
	system = FlexibleSystem(load_system(QE_input_path);
				boundary_conditions=fill(AtomsBase.Periodic(), 3))
	system = attach_psp(system, Dict(atom_symbol => psp_path))
	
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
function run_DFTK_scf(; system, calculator)
	DFTK.reset_timer!(DFTK.timer)
	energy = AtomsCalculators.potential_energy(system, calculator)
	calculator.state.scfres
end

function run_DFTK_vc_relax(system, calculator)
        linesearch =  BackTracking(c_1= 1e-4, ρ_hi= 0.8, ρ_lo= 0.1, order=2, maxstep=Inf)
        solver = OptimizationOptimJL.LBFGS(; linesearch)
        optim_options = (; solver, f_tol=1e-10, g_tol=1e-5, iterations=30,
                         show_trace=true, store_trace = true, allow_f_increases=true)
	DFTK.reset_timer!(DFTK.timer)
        results = minimize_energy!(system, calculator; procedure = "vc_relax", solver, optim_options...)
        (; results, calculator.state.scfres)
end

""" Run a QuantumEspresso computation from a given input file, 
attaching the specified pseudopotential. 
"""
function run_QE(QE_input_path, psp_path; kgrid, Ecut::Unitful.Energy, output_file=nothing)
	Ecut = ustrip(u"Ry", Ecut)
	psp_folder, psp_filename = splitdir(psp_path)
	# Inject correct values for psp and paths. Also inject Ecut and kpoints.
	regex_psp = """s@\\(.*pseudo_dir = \\)\\(.*\\)\$@\\1"$(psp_folder)"/@g"""
	
	regex_kpoints = """s/.*(K_POINTS.*\\n.[^0-9]*)([0-9]*) ([0-9]*) ([0-9]*)(.*)\$/\\1 $(kgrid[1]) $(kgrid[2]) $(kgrid[3])\\5/g"""
	
	regex_Ecut = """s/(.*ecutwfc.[^0-9]*)(.*)/\\1 $Ecut/g"""
	
	input_stream = read(
		pipeline(
		pipeline(
		pipeline(
		`sed -e "$regex_psp" $QE_input_path`,
		`perl -0777 -pe "s/.*(ATOMIC_SPECIES.*\n.[^0-9]*[0-9]*\.[0-9]*)(.*)/\1 $psp_filename/g"`
		),
		`perl -0777 -pe "$regex_kpoints"`
		),
		`perl -0777 -pe "$regex_Ecut"`
		),
		String)
	output_stream = read(pipeline(`echo $input_stream`, `pw.x`), String)
	if !isnothing(output_file)
		run(`echo $output_stream > $output_file`)
	end
	output_stream
end
