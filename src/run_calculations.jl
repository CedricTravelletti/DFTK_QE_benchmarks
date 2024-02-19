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
				   kgrid, Ecut, temperature=1e-4, fft_size_algorithm=:fast)
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

	# Allow for different computations of the FFT size, 
	# in order to diagnose mismatch in number of kpoints
	if fft_size_algorithm == :precise 
		model = model_DFT(system; model_kwargs...)
		
		# This is copied from DFTK/PlaneWaveBasis.jl
		denominators = [denominator(rationalize(sym.w[i]; tol=DFTK.SYMMETRY_TOLERANCE))
                            for sym in model.symmetries for i = 1:3]
		# factors = intersect((2, 3, 4, 6), denominators)
		factors = (1, )
		
		fft_size = compute_fft_size(model, Ecut, kgrid; factors)
	elseif fft_size_algorithm == :fast
		fft_size = nothing
	else
		throw(ErrorException("Unrecognized fft_size algorithm."))
	end
	basis_kwargs = (; kgrid, Ecut, fft_size)
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
		write(output_file, output_stream)
	end
	output_stream
end

""" Run QE vs DFTK comparison on a given system.

Runs both an LDA and a GGA comparison.
User has to provide a QE input file defining the system, as well as pseudopotential files. 
Currently, only single element systems are supported, and the user needs to specify the 
atomic symbol of the element. 
"""
function run_system_comparison(QE_input_path, lda_psp_path, gga_psp_path; atom_symbol, kgrid, Ecut)
	# Name QE output files with as similar name as the input file.
	outfilename_lda = splitdir(QE_input_path)[2][1:end-3] * ".lda.out"
	outfilename_gga = splitdir(QE_input_path)[2][1:end-3] * ".gga.out"
	QE_output_lda = process_QE_output(
				run_QE(QE_input_path, lda_psp_path;
				       kgrid, Ecut,
				       output_file=joinpath(@__DIR__, outfilename_lda)
				       ))
	
	QE_output_gga = process_QE_output(
				run_QE(QE_input_path, gga_psp_path;
				       kgrid, Ecut,
				       output_file=joinpath(@__DIR__, outfilename_gga)
				       ))
	
	# Run with the precise FFT size algorithm to compare number of k-points.
	DFTK_scfres_lda_fft_precise = run_DFTK_scf(;
				build_computation_from_QE(
							  QE_input_path, lda_psp_path;
							  atom_symbol, functional="lda",
							  kgrid, Ecut, fft_size_algorithm=:precise)...)
	DFTK_output_lda_fft_precise = process_DFTK_output(DFTK_scfres_lda_fft_precise)
	
	DFTK_scfres_lda = run_DFTK_scf(;
				build_computation_from_QE(
							  QE_input_path, lda_psp_path;
							  atom_symbol, functional="lda",
							  kgrid, Ecut)...)
	DFTK_output_lda = process_DFTK_output(DFTK_scfres_lda)
	
	
	DFTK_scfres_gga = run_DFTK_scf(;
				build_computation_from_QE(
							  QE_input_path, gga_psp_path;
							  atom_symbol, functional="gga",
							  kgrid, Ecut)...)
	DFTK_output_gga = process_DFTK_output(DFTK_scfres_gga)
	
	
	DFTK_scfres_gga_fft_precise = run_DFTK_scf(;
				build_computation_from_QE(
							  QE_input_path, gga_psp_path;
							  atom_symbol, functional="gga",
							  kgrid, Ecut, fft_size_algorithm=:precise)...)
	DFTK_output_gga_fft_precise = process_DFTK_output(DFTK_scfres_gga_fft_precise)

	(; QE_output_lda, QE_output_gga, DFTK_output_lda, DFTK_output_lda_fft_precise, DFTK_output_gga, DFTK_output_gga_fft_precise)
end
