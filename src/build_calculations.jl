using AtomsBase
using DFTK 


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
