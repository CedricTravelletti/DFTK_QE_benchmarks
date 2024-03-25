using MPI
using TimerOutputs
using DFTK
using AtomsCalculators
using AtomsIO
using AtomsIOPython
using AtomsBase
using Unitful
using UnitfulAtomic


# Helper function to check whether we are on the master process
mpi_master(comm=MPI.COMM_WORLD) = (MPI.Init(); MPI.Comm_rank(comm) == 0)
mpi_nprocs(comm=MPI.COMM_WORLD) = (MPI.Init(); MPI.Comm_size(comm))


""" Run energy computation for using DFTK. 
"""
function run_DFTK_scf(; system, calculator)
	DFTK.reset_timer!(DFTK.timer)
	energy = AtomsCalculators.potential_energy(system, calculator)
	# Dump timings
    	if mpi_master()
    	    println(DFTK.timer)
	    timings = TimerOutputs.todict(DFTK.timer)
        else timings = nothing
	end
	(; calculator.state.scfres, timings)
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
function run_QE(QE_input_path, psp_path; kgrid, Ecut::Unitful.Energy, output_file)
	Ecut = ustrip(u"Ry", Ecut)
	psp_folder, psp_filename = splitdir(psp_path)
	# Inject correct values for psp and paths. Also inject Ecut and kpoints.
	regex_psp = """s@\\(.*pseudo_dir = \\)\\(.*\\)\$@\\1"$(psp_folder)"/@g"""
	regex_outfile = """s@\\(.*pseudo_dir = \\)\\(.*\\)\$@\\1"$(output_file)"/@g"""
	
	regex_kpoints = """s/.*(K_POINTS.*\\n.[^0-9]*)([0-9]*) ([0-9]*) ([0-9]*)(.*)\$/\\1 $(kgrid[1]) $(kgrid[2]) $(kgrid[3])\\5/g"""
	
	regex_Ecut = """s/(.*ecutwfc.[^0-9]*)(.*)/\\1 $Ecut/g"""
	
	
	input_stream = read(
		pipeline(
		pipeline(
		pipeline(
                pipeline(
		`sed -e "$regex_psp" $QE_input_path`,
		`perl -0777 -pe "$regex_outfile"`
		),
		`perl -0777 -pe "s/.*(ATOMIC_SPECIES.*\n.[^0-9]*[0-9]*\.[0-9]*)(.*)/\1 $psp_filename/g"`
		),
		`perl -0777 -pe "$regex_kpoints"`
		),
		`perl -0777 -pe "$regex_Ecut"`
		),
		String)
	output_stream = read(pipeline(`echo $input_stream`, `pw.x`), String)
	output_stream
end

""" Run QE vs DFTK comparison on a given system.

Runs both an LDA and a GGA comparison.
User has to provide a QE input file defining the system, as well as pseudopotential files. 
Currently, only single element systems are supported, and the user needs to specify the 
atomic symbol of the element. 
"""
function run_system_scf_comparison(QE_input_path, lda_psp_path, gga_psp_path, output_folder;
			           atom_symbol, kgrid, Ecut)
	# Name QE output files with as similar name as the input file.
	outfilename_lda = splitdir(QE_input_path)[2][1:end-3] * ".lda.out"
	outfilename_gga = splitdir(QE_input_path)[2][1:end-3] * ".gga.out"
	QE_output_lda = process_QE_output(
				run_QE(QE_input_path, lda_psp_path;
				       kgrid, Ecut,
				       output_file=joinpath(output_folder, outfilename_lda)
				       ))
	
	QE_output_gga = process_QE_output(
				run_QE(QE_input_path, gga_psp_path;
				       kgrid, Ecut,
				       output_file=joinpath(output_folder, outfilename_gga)
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
