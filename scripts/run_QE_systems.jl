using DFTKbenchmarks
using LazyArtifacts
using Unitful
using UnitfulAtomic

Si_QE_input_path = joinpath(@__DIR__, "..", "data/si.scf.in")

# Run both an LDA and a GGA calculation.
Si_lda_psp_path = artifact"pd_nc_sr_lda_standard_0.4.1_upf/Si.upf"
Si_gga_psp_path = artifact"pd_nc_sr_pbe_standard_0.4.1_upf/Si.upf"

QE_output_string = process_QE_output(
			run_QE(Si_QE_input_path, Si_lda_psp_path;
			       k_points=[6, 6, 6], Ecut=30u"Ry"))
QE_output_gga = process_QE_output(
			run_QE(Si_QE_input_path, Si_gga_psp_path;
			       k_points=[6, 6, 6], Ecut=30u"Ry"))
