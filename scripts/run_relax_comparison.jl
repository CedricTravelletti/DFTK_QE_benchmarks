# Compare QE and DFTK on atomic position relaxation (fixed unit cell).
using DFTKbenchmarks
using LazyArtifacts
using Unitful
using UnitfulAtomic


# Silicon
QE_input_path_Si = joinpath(@__DIR__, "..", "data/Si.relax.in")
lda_psp_path_Si = artifact"pd_nc_sr_lda_standard_0.4.1_upf/Si.upf"
gga_psp_path_Si = artifact"pd_nc_sr_pbe_standard_0.4.1_upf/Si.upf"

atom_symbol = :Si
kgrid=[6, 6, 6]
Ecut=50u"Ry"

out_Si = run_system_relax_comparison(QE_input_path_Si, lda_psp_path_Si, gga_psp_path_Si; atom_symbol, kgrid, Ecut)
