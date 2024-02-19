using DFTKbenchmarks
using LazyArtifacts
using Unitful
using UnitfulAtomic


QE_input_path_graphene = joinpath(@__DIR__, "..", "data/graphene.in")
lda_psp_path_graphene = artifact"pd_nc_sr_lda_standard_0.4.1_upf/C.upf"
gga_psp_path_graphene = artifact"pd_nc_sr_pbe_standard_0.4.1_upf/C.upf"

atom_symbol = :C
kgrid=[9, 9, 1]
Ecut=50u"Ry"
	
out_graphene = run_system_comparison(QE_input_path_graphene, lda_psp_path_graphene, gga_psp_path_graphene; atom_symbol, kgrid, Ecut)

# Silicon
QE_input_path_Si = joinpath(@__DIR__, "..", "data/si_alternate.in")
lda_psp_path_Si = artifact"pd_nc_sr_lda_standard_0.4.1_upf/Si.upf"
gga_psp_path_Si = artifact"pd_nc_sr_pbe_standard_0.4.1_upf/Si.upf"

atom_symbol = :Si
kgrid=[6, 6, 6]
Ecut=50u"Ry"

out_Si = run_system_comparison(QE_input_path_Si, lda_psp_path_Si, gga_psp_path_Si; atom_symbol, kgrid, Ecut)
