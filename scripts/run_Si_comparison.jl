using DFTKbenchmarks
using LazyArtifacts
using Unitful
using UnitfulAtomic

Si_QE_input_path = joinpath(@__DIR__, "..", "data/si_alternate.in")

# Run both an LDA and a GGA calculation.
Si_lda_psp_path = artifact"pd_nc_sr_lda_standard_0.4.1_upf/Si.upf"
Si_gga_psp_path = artifact"pd_nc_sr_pbe_standard_0.4.1_upf/Si.upf"

kgrid=[6, 6, 6]
Ecut=50u"Ry"

QE_output_Si_lda = process_QE_output(
			run_QE(Si_QE_input_path, Si_lda_psp_path;
			       kgrid, Ecut,
			       output_file=joinpath(@__DIR__, "tmp/si_scf_lda_QE.out")
			       ))

QE_output_Si_gga = process_QE_output(
			run_QE(Si_QE_input_path, Si_gga_psp_path;
			       kgrid, Ecut,
			       output_file=joinpath(@__DIR__, "tmp/si_scf_gga_QE.out")
			       ))

# Run with the precise FFT size algorithm to compare number of k-points.
DFTK_scfres_Si_lda_fft_precise = run_DFTK_scf(;
			build_computation_from_QE(
						  Si_QE_input_path, Si_lda_psp_path;
						  atom_symbol=:Si, functional="lda",
						  kgrid, Ecut, fft_size_algorithm=:precise)...)
DFTK_output_Si_lda_fft_precise = process_DFTK_output(DFTK_scfres_Si_lda_fft_precise)

DFTK_scfres_Si_lda = run_DFTK_scf(;
			build_computation_from_QE(
						  Si_QE_input_path, Si_lda_psp_path;
						  atom_symbol=:Si, functional="lda",
						  kgrid, Ecut)...)
DFTK_output_Si_lda = process_DFTK_output(DFTK_scfres_Si_lda)


DFTK_scfres_Si_gga = run_DFTK_scf(;
			build_computation_from_QE(
						  Si_QE_input_path, Si_gga_psp_path;
						  atom_symbol=:Si, functional="gga",
						  kgrid, Ecut)...)
DFTK_output_Si_gga = process_DFTK_output(DFTK_scfres_Si_gga)


DFTK_scfres_Si_gga_fft_precise = run_DFTK_scf(;
			build_computation_from_QE(
						  Si_QE_input_path, Si_gga_psp_path;
						  atom_symbol=:Si, functional="gga",
						  kgrid, Ecut, fft_size_algorithm=:precise)...)
DFTK_output_Si_gga_fft_precise = process_DFTK_output(DFTK_scfres_Si_gga_fft_precise)
