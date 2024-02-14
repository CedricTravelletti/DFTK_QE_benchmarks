using DFTKbenchmarks
using LazyArtifacts
using Unitful
using UnitfulAtomic

QE_input_path = joinpath(@__DIR__, "..", "data/graphene.in")

# Run both an LDA and a GGA calculation.
lda_psp_path = artifact"pd_nc_sr_lda_standard_0.4.1_upf/C.upf"
gga_psp_path = artifact"pd_nc_sr_pbe_standard_0.4.1_upf/C.upf"

kgrid=[9, 9, 1]
Ecut=50u"Ry"

QE_output_lda = process_QE_output(
			run_QE(QE_input_path, lda_psp_path;
			       kgrid, Ecut,
			       output_file=joinpath(@__DIR__, "tmp/si_scf_lda_QE.out")
			       ))

QE_output_gga = process_QE_output(
			run_QE(QE_input_path, gga_psp_path;
			       kgrid, Ecut,
			       output_file=joinpath(@__DIR__, "tmp/si_scf_gga_QE.out")
			       ))

# Run with the precise FFT size algorithm to compare number of k-points.
DFTK_scfres_lda_fft_precise = run_DFTK_scf(;
			build_computation_from_QE(
						  QE_input_path, lda_psp_path;
						  atom_symbol=:C, functional="lda",
						  kgrid, Ecut, fft_size_algorithm=:precise)...)
DFTK_output_lda_fft_precise = process_DFTK_output(DFTK_scfres_lda_fft_precise)

DFTK_scfres_lda = run_DFTK_scf(;
			build_computation_from_QE(
						  QE_input_path, lda_psp_path;
						  atom_symbol=:C, functional="lda",
						  kgrid, Ecut)...)
DFTK_output_lda = process_DFTK_output(DFTK_scfres_lda)


DFTK_scfres_gga = run_DFTK_scf(;
			build_computation_from_QE(
						  QE_input_path, gga_psp_path;
						  atom_symbol=:C, functional="gga",
						  kgrid, Ecut)...)
DFTK_output_gga = process_DFTK_output(DFTK_scfres_Si_gga)


DFTK_scfres_gga_fft_precise = run_DFTK_scf(;
			build_computation_from_QE(
						  QE_input_path, gga_psp_path;
						  atom_symbol=:C, functional="gga",
						  kgrid, Ecut, fft_size_algorithm=:precise)...)
DFTK_output_gga_fft_precise = process_DFTK_output(DFTK_scfres_gga_fft_precise)
