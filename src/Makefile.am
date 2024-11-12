lib_LTLIBRARIES= libpsns.la
libpsns_la_SOURCES =  thermal.c diffemissivity.c pair.c photo.c weak.c pair_tbl.c pair_tbl1.c pair_tbl2.c pair_tbl3.c plasma_tbl.c chemical_potential_tbl.c thermal.h thermal_tbl.h weak.h constants.h
#libpsns_CFLAGS=$(OMPFLAGS)

# Can't link Cuba  library into .so forced to use static libs
libpsns_la_LDFLAGS = -static 

EXTRA_DIST=ffn_data.c ffn_table_nue.c ffn_table_nue_A.c ffn_table_nue_B.c ffn_table_nuebar.c ffn_table_nuebar_A.c ffn_table_nuebar_B.c pair_R_table.c pair_Q_table.c pair_table_nue.c pair_table_nuebar.c pair_table_numu.c plasmaL_Q_table.c plasmaL_R_table.c plasmaT_Q_table.c plasmaT_R_table.c chemical_potential_table.c 
